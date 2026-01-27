// ============================================================================
// InterpolationManager.cc - Implementation of solution transfer between meshes
// ============================================================================
// Enhanced version with:
// - Line extraction for pressure and displacement
// - Polynomial fitting using OLS with GMM (no Eigen dependency)
// - Multi-field support (pressure + 3 displacement components)
// - BC callback creation for PLC outer wall
// ============================================================================

#include "../include/InterpolationManager.h"
#include <getfem/getfem_interpolation.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>

// ============================================================================
// Constructor
// ============================================================================
InterpolationManager::InterpolationManager(const GetPot& dataFile,
                                           Bulk* sourceBulk,
                                           Bulk* targetBulk)
    : M_section("russian_doll/interpolation/"),
      M_sourceBulk(sourceBulk),
      M_targetBulk(targetBulk),
      M_z_min(0.0),
      M_z_max(1.0),
      M_matrixBuilt(false)
{
    // Read approach from data file
    std::string approach_str = dataFile((M_section + "approach").c_str(), "line");
    
    if (approach_str == "mesh" || approach_str == "MESH") {
        M_approach = Approach::MESH_INTERPOLATION;
    } else {
        M_approach = Approach::LINE_INTERPOLATION;
        readLineProfilesFromFile(dataFile);
    }
    
    std::cout << "[InterpolationManager] Created with approach: " 
              << (M_approach == Approach::LINE_INTERPOLATION ? "LINE" : "MESH") 
              << std::endl;
}

// ============================================================================
// Initialize
// ============================================================================
void InterpolationManager::initialize() {
    if (M_approach == Approach::LINE_INTERPOLATION) {
        if (M_lineProfiles.empty()) {
            std::cerr << "[InterpolationManager] Warning: No line profiles defined!" << std::endl;
        } else {
            std::cout << "[InterpolationManager] Initialized with " 
                      << M_lineProfiles.size() << " line profiles" << std::endl;
            
            // Print line profile info
            for (const auto& profile : M_lineProfiles) {
                std::cout << "  Line '" << profile.name << "':" << std::endl;
                std::cout << "    Start: (" << profile.start_point[0] << ", " 
                          << profile.start_point[1] << ", " << profile.start_point[2] << ")" << std::endl;
                std::cout << "    End: (" << profile.end_point[0] << ", " 
                          << profile.end_point[1] << ", " << profile.end_point[2] << ")" << std::endl;
                std::cout << "    Samples: " << profile.num_samples << std::endl;
                std::cout << "    Polynomial order: " << profile.polynomial_order << std::endl;
            }
        }
    } else {
        std::cout << "[InterpolationManager] MESH approach - will build matrix on first use" << std::endl;
    }
}

// ============================================================================
// Read Line Profiles from Data File
// ============================================================================
void InterpolationManager::readLineProfilesFromFile(const GetPot& dataFile) {
    size_type num_lines = dataFile((M_section + "num_lines").c_str(), 1);
    
    for (size_type i = 0; i < num_lines; ++i) {
        std::string prefix = M_section + "line" + std::to_string(i) + "/";
        
        LineProfile profile;
        profile.name = dataFile((prefix + "name").c_str(), 
                                ("line" + std::to_string(i)).c_str());
        
        // Read start point
        profile.start_point.resize(3);
        profile.start_point[0] = dataFile((prefix + "start_x").c_str(), 0.0);
        profile.start_point[1] = dataFile((prefix + "start_y").c_str(), 0.0);
        profile.start_point[2] = dataFile((prefix + "start_z").c_str(), 0.0);
        
        // Read end point
        profile.end_point.resize(3);
        profile.end_point[0] = dataFile((prefix + "end_x").c_str(), 0.0);
        profile.end_point[1] = dataFile((prefix + "end_y").c_str(), 0.0);
        profile.end_point[2] = dataFile((prefix + "end_z").c_str(), 1.0);
        
        // Sampling parameters
        profile.num_samples = dataFile((prefix + "num_samples").c_str(), 100);
        profile.polynomial_order = dataFile((prefix + "polynomial_order").c_str(), 3);
        
        M_lineProfiles.push_back(profile);
        
        std::cout << "[InterpolationManager] Loaded line profile '" << profile.name 
                  << "' from data file" << std::endl;
    }
}

// ============================================================================
// Public Methods
// ============================================================================
void InterpolationManager::addLineProfile(const LineProfile& profile) {
    M_lineProfiles.push_back(profile);
    std::cout << "[InterpolationManager] Added line profile '" << profile.name << "'" << std::endl;
}

// ============================================================================
// Helper: Interpolate Point Along Line
// ============================================================================
bgeot::base_node InterpolationManager::interpolatePoint(const LineProfile& profile, 
                                                        scalar_type t) {
    bgeot::base_node pt(3);
    for (int d = 0; d < 3; ++d) {
        pt[d] = (1.0 - t) * profile.start_point[d] + t * profile.end_point[d];
    }
    return pt;
}

// ============================================================================
// Helper: Compute Z-coordinates from sample points
// ============================================================================
void InterpolationManager::computeZCoordinates(const std::vector<bgeot::base_node>& sample_points,
                                               std::vector<scalar_type>& z_coords) {
    z_coords.resize(sample_points.size());
    for (size_type i = 0; i < sample_points.size(); ++i) {
        z_coords[i] = sample_points[i][2];  // Assuming z is the 3rd coordinate
    }
}

// ============================================================================
// Extract Scalar Solution Along a Line
// ============================================================================
void InterpolationManager::extractAlongLine(const scalarVectorPtr_Type& solution,
                                            const getfem::mesh_fem& mf_source,
                                            const LineProfile& profile,
                                            std::vector<scalar_type>& arc_coords,
                                            std::vector<scalar_type>& values) {
    arc_coords.clear();
    values.clear();
    arc_coords.reserve(profile.num_samples);
    values.reserve(profile.num_samples);
    
    // Initialize z-range
    M_z_min = std::numeric_limits<scalar_type>::max();
    M_z_max = std::numeric_limits<scalar_type>::lowest();
    
    // Compute total arc length and direction
    bgeot::base_node diff = profile.end_point - profile.start_point;
    scalar_type total_length = gmm::vect_norm2(diff);
    
    if (total_length < 1e-15) {
        std::cerr << "[InterpolationManager] Error: Zero-length line profile!" << std::endl;
        return;
    }
    
    // Store sample points with their parametric coordinate t
    std::vector<bgeot::base_node> sample_points(profile.num_samples);
    std::vector<scalar_type> sample_t(profile.num_samples);
    
    // Create 1D mesh with points along the line
    getfem::mesh mesh_1d;
    std::vector<size_type> pt_ids;
    pt_ids.reserve(profile.num_samples);
    
    for (size_type i = 0; i < profile.num_samples; ++i) {
        scalar_type t = static_cast<scalar_type>(i) / (profile.num_samples - 1);
        sample_t[i] = t;
        sample_points[i] = interpolatePoint(profile, t);
        pt_ids.push_back(mesh_1d.add_point(sample_points[i]));
        
        // Update z-range (assuming z is the 3rd coordinate)
        scalar_type z = sample_points[i][2];
        if (z < M_z_min) M_z_min = z;
        if (z > M_z_max) M_z_max = z;
    }
    
    // Add 1D segments connecting consecutive points
    for (size_type i = 0; i < profile.num_samples - 1; ++i) {
        std::vector<size_type> ind = {pt_ids[i], pt_ids[i+1]};
        mesh_1d.add_convex(bgeot::simplex_geotrans(1, 1), ind.begin());
    }
    
    // Create mesh_fem on 1D mesh (P1 Lagrange - DOFs at vertices)
    getfem::mesh_fem mf_1d(mesh_1d, 1);  // qdim = 1 (scalar)
    mf_1d.set_classical_finite_element(1);  // P1 elements
    
    size_type nb_dof_1d = mf_1d.nb_dof();
    
    // Allocate interpolation target vector
    std::vector<scalar_type> interp_values(nb_dof_1d, 0.0);
    
    try {
        // Perform 3D -> 1D interpolation using GetFEM
        getfem::interpolation(mf_source, mf_1d, *solution, interp_values);
    } catch (std::exception& e) {
        std::cerr << "[InterpolationManager] Interpolation failed: " << e.what() << std::endl;
        return;
    }
    
    // Map DOFs back to sample points using DOF coordinates
    values.resize(profile.num_samples);
    arc_coords.resize(profile.num_samples);
    
    // For each sample point, find the corresponding DOF
    for (size_type i = 0; i < profile.num_samples; ++i) {
        arc_coords[i] = sample_t[i] * total_length;
        
        // Find DOF closest to this sample point
        scalar_type min_dist = std::numeric_limits<scalar_type>::max();
        size_type best_dof = 0;
        
        for (size_type d = 0; d < nb_dof_1d; ++d) {
            bgeot::base_node dof_pt = mf_1d.point_of_basic_dof(d);
            scalar_type dist = gmm::vect_dist2(sample_points[i], dof_pt);
            if (dist < min_dist) {
                min_dist = dist;
                best_dof = d;
            }
        }
        
        values[i] = interp_values[best_dof];
    }
    #ifdef VERBOSE
    std::cout << "[InterpolationManager] Extracted " << values.size() 
              << " samples along line '" << profile.name 
              << "' (length = " << total_length 
              << ", z-range = [" << M_z_min << ", " << M_z_max << "])" << std::endl;
    #endif
}

// ============================================================================
// Extract Vector Solution (Displacement) Along a Line
// ============================================================================
void InterpolationManager::extractVectorAlongLine(const scalarVectorPtr_Type& solution,
                                                  const getfem::mesh_fem& mf_source,
                                                  const LineProfile& profile,
                                                  std::vector<scalar_type>& arc_coords,
                                                  std::vector<scalar_type>& values_x,
                                                  std::vector<scalar_type>& values_y,
                                                  std::vector<scalar_type>& values_z) {
    arc_coords.clear();
    values_x.clear();
    values_y.clear();
    values_z.clear();
    
    arc_coords.reserve(profile.num_samples);
    values_x.reserve(profile.num_samples);
    values_y.reserve(profile.num_samples);
    values_z.reserve(profile.num_samples);
    
    // Compute total arc length
    bgeot::base_node diff = profile.end_point - profile.start_point;
    scalar_type total_length = gmm::vect_norm2(diff);
    
    if (total_length < 1e-15) {
        std::cerr << "[InterpolationManager] Error: Zero-length line profile!" << std::endl;
        return;
    }
    
    // Store sample points
    std::vector<bgeot::base_node> sample_points(profile.num_samples);
    std::vector<scalar_type> sample_t(profile.num_samples);
    
    // Create 1D mesh
    getfem::mesh mesh_1d;
    std::vector<size_type> pt_ids;
    pt_ids.reserve(profile.num_samples);
    
    for (size_type i = 0; i < profile.num_samples; ++i) {
        scalar_type t = static_cast<scalar_type>(i) / (profile.num_samples - 1);
        sample_t[i] = t;
        sample_points[i] = interpolatePoint(profile, t);
        pt_ids.push_back(mesh_1d.add_point(sample_points[i]));
    }
    
    // Add 1D segments
    for (size_type i = 0; i < profile.num_samples - 1; ++i) {
        std::vector<size_type> ind = {pt_ids[i], pt_ids[i+1]};
        mesh_1d.add_convex(bgeot::simplex_geotrans(1, 1), ind.begin());
    }
    
    // Create mesh_fem on 1D mesh with qdim = 3 (vector field)
    getfem::mesh_fem mf_1d(mesh_1d, 3);  // qdim = 3 for vector
    mf_1d.set_classical_finite_element(1);
    
    size_type nb_dof_1d = mf_1d.nb_dof();  // This will be 3 * num_points
    
    // Allocate interpolation target vector
    std::vector<scalar_type> interp_values(nb_dof_1d, 0.0);
    
    try {
        getfem::interpolation(mf_source, mf_1d, *solution, interp_values);
    } catch (std::exception& e) {
        std::cerr << "[InterpolationManager] Vector interpolation failed: " << e.what() << std::endl;
        return;
    }
    
    // Extract components
    values_x.resize(profile.num_samples);
    values_y.resize(profile.num_samples);
    values_z.resize(profile.num_samples);
    arc_coords.resize(profile.num_samples);
    
    size_type nb_pts = profile.num_samples;
    
    for (size_type i = 0; i < nb_pts; ++i) {
        arc_coords[i] = sample_t[i] * total_length;
        
        // Find DOF indices for this point
        // For vector FEM with qdim=3, DOFs are typically ordered by component
        // x-components at indices 0..N-1, y at N..2N-1, z at 2N..3N-1
        // Or they might be interleaved: [u0x, u0y, u0z, u1x, u1y, u1z, ...]
        // GetFEM typically uses interleaved format
        
        // Find closest point's index
        scalar_type min_dist = std::numeric_limits<scalar_type>::max();
        size_type best_pt = 0;
        
        for (size_type d = 0; d < nb_pts; ++d) {
            bgeot::base_node dof_pt = mf_1d.point_of_basic_dof(3*d);  // x-component DOF
            scalar_type dist = gmm::vect_dist2(sample_points[i], dof_pt);
            if (dist < min_dist) {
                min_dist = dist;
                best_pt = d;
            }
        }
        
        // Extract interleaved components
        values_x[i] = interp_values[3*best_pt];
        values_y[i] = interp_values[3*best_pt + 1];
        values_z[i] = interp_values[3*best_pt + 2];
    }
    
    std::cout << "[InterpolationManager] Extracted vector field (" << profile.num_samples 
              << " samples) along line '" << profile.name << "'" << std::endl;
}

// ============================================================================
// GMM Helper: Solve Normal Equations
// ============================================================================
void InterpolationManager::solveNormalEquations(const gmm::dense_matrix<scalar_type>& VtV,
                                                const std::vector<scalar_type>& Vty,
                                                std::vector<scalar_type>& coeffs) {
    size_type n = gmm::mat_nrows(VtV);
    
    // Copy VtV because LU modifies it
    gmm::dense_matrix<scalar_type> A(n, n);
    gmm::copy(VtV, A);
    
    // Add small regularization for numerical stability
    for (size_type i = 0; i < n; ++i) {
        A(i, i) += 1.0e-12;
    }
    
    // LU factorization
    std::vector<size_type> ipvt(n);
    gmm::lu_factor(A, ipvt);
    
    // Solve
    coeffs.resize(n);
    gmm::copy(Vty, coeffs);
    gmm::lu_solve(A, ipvt, coeffs, Vty);
}

// ============================================================================
// Fit Polynomial Using OLS (Ordinary Least Squares)
// ============================================================================
PolynomialFit InterpolationManager::fitPolynomial(const std::vector<scalar_type>& arc_coords,
                                                  const std::vector<scalar_type>& values,
                                                  size_type order,
                                                  const std::string& field_name) {
    PolynomialFit result;
    result.field_name = field_name;
    result.z_min = M_z_min;
    result.z_max = M_z_max;
    
    size_type n = values.size();
    size_type m = order + 1;  // Number of coefficients
    
    if (n < m) {
        std::cerr << "[InterpolationManager] Error: Not enough points (" << n 
                  << ") for polynomial order " << order << std::endl;
        result.coefficients.resize(m, 0.0);
        return result;
    }
    
    // Normalize arc coordinates to [0, 1]
    scalar_type max_arc = arc_coords.back();
    result.arc_length = max_arc;
    if (max_arc < 1.0e-15) max_arc = 1.0;
    
    // Build Vandermonde matrix V (n x m) using GMM dense matrix
    gmm::dense_matrix<scalar_type> V(n, m);
    gmm::clear(V);
    
    for (size_type i = 0; i < n; ++i) {
        scalar_type t = arc_coords[i] / max_arc;  // Normalized coordinate
        scalar_type t_power = 1.0;
        for (size_type j = 0; j < m; ++j) {
            V(i, j) = t_power;
            t_power *= t;
        }
    }
    
    // Compute V^T * V (m x m)
    gmm::dense_matrix<scalar_type> VtV(m, m);
    gmm::clear(VtV);
    gmm::mult(gmm::transposed(V), V, VtV);
    
    // Compute V^T * y (m x 1)
    std::vector<scalar_type> Vty(m, 0.0);
    gmm::mult(gmm::transposed(V), values, Vty);
    
    // Solve normal equations: (V^T V) * coeffs = V^T * y
    solveNormalEquations(VtV, Vty, result.coefficients);
    
    // Compute R² goodness of fit
    std::vector<scalar_type> y_pred(n, 0.0);
    gmm::mult(V, result.coefficients, y_pred);
    
    // ss_res = sum((y - y_pred)^2)
    scalar_type ss_res = 0.0;
    for (size_type i = 0; i < n; ++i) {
        scalar_type residual = values[i] - y_pred[i];
        ss_res += residual * residual;
    }
    
    // mean_y
    scalar_type mean_y = 0.0;
    for (size_type i = 0; i < n; ++i) {
        mean_y += values[i];
    }
    mean_y /= n;
    
    // ss_tot = sum((y - mean_y)^2)
    scalar_type ss_tot = 0.0;
    for (size_type i = 0; i < n; ++i) {
        scalar_type diff = values[i] - mean_y;
        ss_tot += diff * diff;
    }
    
    result.r_squared = (ss_tot > 1.0e-15) ? (1.0 - ss_res / ss_tot) : 1.0;
    #ifdef VERBOSE
    std::cout << "[InterpolationManager] Polynomial fit for '" << field_name 
              << "' (order " << order << "): R² = " << result.r_squared << std::endl;
    std::cout << "  Coefficients: [";

    for (size_type i = 0; i < result.coefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << result.coefficients[i];
    }
    std::cout << "]" << std::endl;
    #endif
    std::cout << "[Interpolation Manager] Polynomial fit completed for " << field_name << std::endl;    
    return result;
}

// ============================================================================
// Fit Polynomial with Specified Z-Range
// ============================================================================
PolynomialFit InterpolationManager::fitPolynomialWithZRange(
    const std::vector<scalar_type>& z_coords,
    const std::vector<scalar_type>& values,
    size_type order,
    scalar_type z_min,
    scalar_type z_max,
    const std::string& field_name) {
    
    PolynomialFit result;
    result.field_name = field_name;
    result.z_min = z_min;
    result.z_max = z_max;
    
    size_type n = values.size();
    size_type m = order + 1;
    
    if (n < m) {
        std::cerr << "[InterpolationManager] Error: Not enough points for polynomial order " 
                  << order << std::endl;
        result.coefficients.resize(m, 0.0);
        return result;
    }
    
    // Build Vandermonde matrix with normalized z
    gmm::dense_matrix<scalar_type> V(n, m);
    gmm::clear(V);
    
    scalar_type dz = z_max - z_min;
    if (std::abs(dz) < 1e-15) dz = 1.0;
    
    for (size_type i = 0; i < n; ++i) {
        scalar_type t = (z_coords[i] - z_min) / dz;  // Normalize to [0,1]
        t = std::max(0.0, std::min(1.0, t));
        
        scalar_type t_power = 1.0;
        for (size_type j = 0; j < m; ++j) {
            V(i, j) = t_power;
            t_power *= t;
        }
    }
    
    // Compute V^T * V
    gmm::dense_matrix<scalar_type> VtV(m, m);
    gmm::clear(VtV);
    gmm::mult(gmm::transposed(V), V, VtV);
    
    // Compute V^T * y
    std::vector<scalar_type> Vty(m, 0.0);
    gmm::mult(gmm::transposed(V), values, Vty);
    
    // Solve
    solveNormalEquations(VtV, Vty, result.coefficients);
    
    // Compute R²
    std::vector<scalar_type> y_pred(n, 0.0);
    gmm::mult(V, result.coefficients, y_pred);
    
    scalar_type ss_res = 0.0;
    scalar_type mean_y = 0.0;
    for (size_type i = 0; i < n; ++i) {
        scalar_type residual = values[i] - y_pred[i];
        ss_res += residual * residual;
        mean_y += values[i];
    }
    mean_y /= n;
    
    scalar_type ss_tot = 0.0;
    for (size_type i = 0; i < n; ++i) {
        scalar_type diff = values[i] - mean_y;
        ss_tot += diff * diff;
    }
    
    result.r_squared = (ss_tot > 1.0e-15) ? (1.0 - ss_res / ss_tot) : 1.0;
    
    return result;
}

// ============================================================================
// Interpolate Pressure
// ============================================================================
PolynomialFit InterpolationManager::interpolatePressure(
    const scalarVectorPtr_Type& pressure_solution,
    const getfem::mesh_fem& mf_pressure) {
    
    if (M_lineProfiles.empty()) {
        std::cerr << "[InterpolationManager] Error: No line profiles defined!" << std::endl;
        return PolynomialFit();
    }
    
    const LineProfile& profile = M_lineProfiles[0];
    
    // Extract pressure along line
    std::vector<scalar_type> arc_coords, values;
    extractAlongLine(pressure_solution, mf_pressure, profile, arc_coords, values);
    
    if (values.empty()) {
        std::cerr << "[InterpolationManager] Error: No pressure values extracted!" << std::endl;
        return PolynomialFit();
    }
    
    // Fit polynomial
    PolynomialFit fit = fitPolynomial(arc_coords, values, profile.polynomial_order, "pressure");
    fit.z_min = M_z_min;
    fit.z_max = M_z_max;
    
    // Store in field polynomials
    M_fieldPolynomials.pressure = fit;
    M_fieldPolynomials.z_min = M_z_min;
    M_fieldPolynomials.z_max = M_z_max;
    
    // Also add to legacy list
    M_polynomialFits.clear();
    M_polynomialFits.push_back(fit);
    
    return fit;
}

// ============================================================================
// Interpolate Displacement
// ============================================================================
void InterpolationManager::interpolateDisplacement(
    const scalarVectorPtr_Type& displacement_solution,
    const getfem::mesh_fem& mf_displacement,
    PolynomialFit& fit_x,
    PolynomialFit& fit_y,
    PolynomialFit& fit_z) {
    
    if (M_lineProfiles.empty()) {
        std::cerr << "[InterpolationManager] Error: No line profiles defined!" << std::endl;
        return;
    }
    
    const LineProfile& profile = M_lineProfiles[0];
    
    // Extract displacement along line
    std::vector<scalar_type> arc_coords;
    std::vector<scalar_type> values_x, values_y, values_z;
    
    extractVectorAlongLine(displacement_solution, mf_displacement, profile,
                           arc_coords, values_x, values_y, values_z);
    
    if (values_x.empty()) {
        std::cerr << "[InterpolationManager] Error: No displacement values extracted!" << std::endl;
        return;
    }
    
    // Fit polynomials for each component
    fit_x = fitPolynomial(arc_coords, values_x, profile.polynomial_order, "displacement_x");
    fit_y = fitPolynomial(arc_coords, values_y, profile.polynomial_order, "displacement_y");
    fit_z = fitPolynomial(arc_coords, values_z, profile.polynomial_order, "displacement_z");
    
    // Set z-range
    fit_x.z_min = M_z_min; fit_x.z_max = M_z_max;
    fit_y.z_min = M_z_min; fit_y.z_max = M_z_max;
    fit_z.z_min = M_z_min; fit_z.z_max = M_z_max;
    
    // Store in field polynomials
    M_fieldPolynomials.displacement_x = fit_x;
    M_fieldPolynomials.displacement_y = fit_y;
    M_fieldPolynomials.displacement_z = fit_z;
}

// ============================================================================
// Interpolate All Fields
// ============================================================================
FieldPolynomials InterpolationManager::interpolateAllFields(
    const scalarVectorPtr_Type& pressure_solution,
    const getfem::mesh_fem& mf_pressure,
    const scalarVectorPtr_Type& displacement_solution,
    const getfem::mesh_fem& mf_displacement) {
    
    std::cout << "[InterpolationManager] Interpolating all fields..." << std::endl;
    
    // Interpolate pressure
    M_fieldPolynomials.pressure = interpolatePressure(pressure_solution, mf_pressure);
    
    // Interpolate displacement
    interpolateDisplacement(displacement_solution, mf_displacement,
                            M_fieldPolynomials.displacement_x,
                            M_fieldPolynomials.displacement_y,
                            M_fieldPolynomials.displacement_z);
    
    M_fieldPolynomials.z_min = M_z_min;
    M_fieldPolynomials.z_max = M_z_max;
    
    std::cout << "[InterpolationManager] All fields interpolated successfully" << std::endl;
    std::cout << "  Z-range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    
    return M_fieldPolynomials;
}

// ============================================================================
// Apply Polynomial as BC on Target Mesh
// ============================================================================
void InterpolationManager::applyPolynomialBC(const PolynomialFit& fit,
                                             const LineProfile& profile,
                                             const getfem::mesh_fem& mf_target,
                                             scalarVectorPtr_Type& bc_values) {
    size_type nb_dof = mf_target.nb_dof();
    bc_values.reset(new scalarVector_Type(nb_dof, 0.0));
    
    // For each DOF, evaluate polynomial based on its z-coordinate
    for (size_type i = 0; i < nb_dof; ++i) {
        bgeot::base_node dof_pt = mf_target.point_of_basic_dof(i);
        scalar_type z = dof_pt[2];
        (*bc_values)[i] = fit.evaluateAtZ(z);
    }
    
    std::cout << "[InterpolationManager] Applied polynomial BC to " << nb_dof 
              << " DOFs" << std::endl;
}

// ============================================================================
// Create BC Callback
// ============================================================================
std::function<scalar_type(scalar_type)> 
InterpolationManager::createBCCallback(const PolynomialFit& fit) {
    // Capture fit by copy
    return [fit](scalar_type z) -> scalar_type {
        return fit.evaluateAtZ(z);
    };
}

// ============================================================================
// Legacy: Line Interpolation Pipeline
// ============================================================================
void InterpolationManager::interpolateViaLines(const scalarVectorPtr_Type& source_solution,
                                               const getfem::mesh_fem& mf_source,
                                               scalarVectorPtr_Type& target_solution,
                                               const getfem::mesh_fem& mf_target) {
    if (M_lineProfiles.empty()) {
        std::cerr << "[InterpolationManager] Error: No line profiles configured!" << std::endl;
        return;
    }
    
    M_polynomialFits.clear();
    
    const LineProfile& profile = M_lineProfiles[0];
    
    // Step 1: Extract values along line
    std::vector<scalar_type> arc_coords, values;
    extractAlongLine(source_solution, mf_source, profile, arc_coords, values);
    
    if (values.empty()) {
        std::cerr << "[InterpolationManager] Error: No values extracted!" << std::endl;
        return;
    }
    
    // Step 2: Fit polynomial
    PolynomialFit fit = fitPolynomial(arc_coords, values, profile.polynomial_order);
    fit.t_min = 0.0;
    fit.t_max = 1.0;
    M_polynomialFits.push_back(fit);
    
    // Step 3: Apply to target mesh
    applyPolynomialBC(fit, profile, mf_target, target_solution);
}

// ============================================================================
// Legacy: Mesh Interpolation
// ============================================================================
void InterpolationManager::interpolateViaMesh(const scalarVectorPtr_Type& source_solution,
                                              const getfem::mesh_fem& mf_source,
                                              scalarVectorPtr_Type& target_solution,
                                              const getfem::mesh_fem& mf_target) {
    size_type nb_dof_target = mf_target.nb_dof();
    target_solution.reset(new scalarVector_Type(nb_dof_target, 0.0));
    
    try {
        getfem::interpolation(mf_source, mf_target, *source_solution, *target_solution);
        std::cout << "[InterpolationManager] Mesh interpolation: " 
                  << mf_source.nb_dof() << " DOFs -> " 
                  << nb_dof_target << " DOFs" << std::endl;
    } catch (std::exception& e) {
        std::cerr << "[InterpolationManager] Mesh interpolation failed: " << e.what() << std::endl;
    }
}

void InterpolationManager::buildInterpolationMatrix(const getfem::mesh_fem& mf_source,
                                                    const getfem::mesh_fem& mf_target) {
    size_type nb_dof_source = mf_source.nb_dof();
    size_type nb_dof_target = mf_target.nb_dof();
    
    gmm::resize(M_interpMatrix, nb_dof_target, nb_dof_source);
    gmm::clear(M_interpMatrix);
    
    try {
        getfem::interpolation(mf_source, mf_target, M_interpMatrix);
        M_matrixBuilt = true;
        std::cout << "[InterpolationManager] Built interpolation matrix: " 
                  << nb_dof_target << " x " << nb_dof_source << std::endl;
    } catch (std::exception& e) {
        std::cerr << "[InterpolationManager] Failed to build interpolation matrix: " 
                  << e.what() << std::endl;
    }
}

void InterpolationManager::interpolateViaMatrix(const scalarVectorPtr_Type& source_solution,
                                                scalarVectorPtr_Type& target_solution) {
    if (!M_matrixBuilt) {
        std::cerr << "[InterpolationManager] Error: Interpolation matrix not built!" << std::endl;
        return;
    }
    
    size_type nb_dof_target = gmm::mat_nrows(M_interpMatrix);
    target_solution.reset(new scalarVector_Type(nb_dof_target));
    
    gmm::mult(M_interpMatrix, *source_solution, *target_solution);
}

void InterpolationManager::interpolate(const scalarVectorPtr_Type& source_solution,
                                       const getfem::mesh_fem& mf_source,
                                       scalarVectorPtr_Type& target_solution,
                                       const getfem::mesh_fem& mf_target) {
    if (M_approach == Approach::LINE_INTERPOLATION) {
        interpolateViaLines(source_solution, mf_source, target_solution, mf_target);
    } else {
        if (!M_matrixBuilt) {
            buildInterpolationMatrix(mf_source, mf_target);
        }
        interpolateViaMatrix(source_solution, target_solution);
    }
}

// ============================================================================
// Coefficient Access Methods
// ============================================================================
const PolynomialFit& InterpolationManager::getLastPressureFit() const {
    return M_fieldPolynomials.pressure;
}

const PolynomialFit& InterpolationManager::getLastPolynomialFit() const {
    if (M_polynomialFits.empty()) {
        static PolynomialFit empty;
        return empty;
    }
    return M_polynomialFits.back();
}

const std::vector<scalar_type>& InterpolationManager::getPolynomialCoefficients() const {
    return M_fieldPolynomials.pressure.coefficients;
}

const std::vector<scalar_type>& InterpolationManager::getPressureCoefficients() const {
    return M_fieldPolynomials.pressure.coefficients;
}

const std::vector<scalar_type>& InterpolationManager::getDisplacementXCoefficients() const {
    return M_fieldPolynomials.displacement_x.coefficients;
}

const std::vector<scalar_type>& InterpolationManager::getDisplacementYCoefficients() const {
    return M_fieldPolynomials.displacement_y.coefficients;
}

const std::vector<scalar_type>& InterpolationManager::getDisplacementZCoefficients() const {
    return M_fieldPolynomials.displacement_z.coefficients;
}

void InterpolationManager::getZRange(scalar_type& z_min, scalar_type& z_max) const {
    z_min = M_z_min;
    z_max = M_z_max;
}

scalar_type InterpolationManager::evaluatePressureAtZ(scalar_type z) const {
    return M_fieldPolynomials.pressure.evaluateAtZ(z);
}

scalar_type InterpolationManager::evaluateAtZ(scalar_type z) const {
    return evaluatePressureAtZ(z);
}

void InterpolationManager::evaluateDisplacementAtZ(scalar_type z, 
                                                   scalar_type& ux, 
                                                   scalar_type& uy, 
                                                   scalar_type& uz) const {
    ux = M_fieldPolynomials.displacement_x.evaluateAtZ(z);
    uy = M_fieldPolynomials.displacement_y.evaluateAtZ(z);
    uz = M_fieldPolynomials.displacement_z.evaluateAtZ(z);
}