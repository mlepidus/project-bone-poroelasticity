// ============================================================================
// InterpolationManager.cpp - Implementation of solution transfer between meshes
// ============================================================================

#include "../include/InterpolationManager.h"
#include <getfem/getfem_interpolation.h>
#include <iostream>
#include <cmath>
#include <limits>

// ============================================================================
// Constructor
// ============================================================================
InterpolationManager::InterpolationManager(const GetPot& dataFile,
                                           Bulk* sourceBulk,
                                           Bulk* targetBulk)
    :  M_section("russian_doll/interpolation/"),
      M_sourceBulk(sourceBulk),
      M_targetBulk(targetBulk),
      M_z_min(0.0),
      M_z_max(1.0),
      M_matrixBuilt(false)
{
    // Read approach from data file
    std::string approach_str = dataFile((M_section + "approach").c_str(), "mesh");
    
    if (approach_str == "line" || approach_str == "LINE") {
        M_approach = Approach::LINE_INTERPOLATION;
        readLineProfilesFromFile(dataFile);
    } else {
        M_approach = Approach::MESH_INTERPOLATION;
    }
    
    std::cout << "[InterpolationManager] Using approach: " 
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
        }
        std::cout << "[InterpolationManager] Initialized with " 
                  << M_lineProfiles.size() << " line profiles" << std::endl;
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
        profile.end_point[0] = dataFile((prefix + "end_x").c_str(), 1.0);
        profile.end_point[1] = dataFile((prefix + "end_y").c_str(), 0.0);
        profile.end_point[2] = dataFile((prefix + "end_z").c_str(), 0.0);
        
        // Sampling parameters
        profile.num_samples = dataFile((prefix + "num_samples").c_str(), 100);
        profile.polynomial_order = dataFile((prefix + "polynomial_order").c_str(), 3);
        
        M_lineProfiles.push_back(profile);
        
        std::cout << "[InterpolationManager] Added line '" << profile.name 
                  << "' from (" << profile.start_point[0] << "," 
                  << profile.start_point[1] << "," << profile.start_point[2] 
                  << ") to (" << profile.end_point[0] << "," 
                  << profile.end_point[1] << "," << profile.end_point[2] << ")" 
                  << std::endl;
    }
}

// ============================================================================
// Public Methods
// ============================================================================
void InterpolationManager::addLineProfile(const LineProfile& profile) {
    M_lineProfiles.push_back(profile);
}


const PolynomialFit& InterpolationManager::getLastPolynomialFit() const {
    if (M_polynomialFits.empty()) {
        throw std::runtime_error("No polynomial fits available");
    }
    return M_polynomialFits.back();
}

const std::vector<scalar_type>& InterpolationManager::getPolynomialCoefficients() const {
    if (M_polynomialFits.empty()) {
        throw std::runtime_error("No polynomial fits available");
    }
    return M_polynomialFits.back().coefficients;
}

void InterpolationManager::getZRange(scalar_type& z_min, scalar_type& z_max) const {
    z_min = M_z_min;
    z_max = M_z_max;
}

scalar_type InterpolationManager::evaluateAtZ(scalar_type z) const {
    if (M_polynomialFits.empty()) {
        throw std::runtime_error("No polynomial fits available");
    }
    return M_polynomialFits.back().evaluateAtZ(z, M_z_min, M_z_max);
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
// Approach 1: Extract Solution Along a Line (FIXED VERSION)
// ============================================================================
// ============================================================================
// Approach 1: Extract Solution Along a Line (FIXED VERSION)
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
    
    // -------------------------------------------------------------------------
    // Method: Create a temporary 1D mesh along the line and use GetFEM interp
    // -------------------------------------------------------------------------
    
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
    
    // Verify DOF count matches sample count (should be true for P1)
    size_type nb_dof_1d = mf_1d.nb_dof();
    if (nb_dof_1d != profile.num_samples) {
        std::cerr << "[InterpolationManager] Warning: DOF count (" << nb_dof_1d 
                  << ") != sample count (" << profile.num_samples << ")" << std::endl;
    }
    
    /*
    // Check if points are inside the source mesh (diagnostic)
    std::cout << "[InterpolationManager] Checking if line is inside PV mesh..." << std::endl;
    size_type points_found = 0;
    
    // Get the mesh from the mesh_fem
    const getfem::mesh& source_mesh = mf_source.linked_mesh();
    
    for (size_type i = 0; i < profile.num_samples; ++i) {
        bgeot::base_node pt = sample_points[i];
        
        // Method 1: Using convex_of_point - returns -1 if not found
        size_type cv = source_mesh.convex_of_point(pt);
        bool found = (cv != size_type(-1));
        
        // Alternative Method 2: If convex_of_point doesn't work well, you can try:
        // bool found = false;
        // for (getfem::mr_visitor cv(source_mesh.convex_index()); !cv.finished(); ++cv) {
        //     if (source_mesh.convex_contains_point(cv.cv(), pt)) {
        //         found = true;
        //         break;
        //     }
        // }
        
        if (found) {
            points_found++;
        } else if (i % 10 == 0) {  // Print every 10th point
            std::cerr << "[InterpolationManager] WARNING: Point " << i << " at (" 
                      << pt[0] << "," << pt[1] << "," << pt[2] 
                      << ") not in PV mesh!" << std::endl;
        }
    }
    
    std::cout << "[InterpolationManager] " << points_found << " / " << profile.num_samples 
              << " points found in PV mesh" << std::endl;
    
    if (points_found == 0) {
        std::cerr << "[InterpolationManager] ERROR: No points found in PV mesh!" << std::endl;
        return;
    }
    */
    // Allocate interpolation target vector
    std::vector<scalar_type> interp_values(nb_dof_1d, 0.0);
    
    try {
        // Perform 3D -> 1D interpolation using GetFEM
        getfem::interpolation(mf_source, mf_1d, *solution, interp_values);
    } catch (std::exception& e) {
        std::cerr << "[InterpolationManager] Interpolation failed: " << e.what() << std::endl;
        return;
    }
    
    // -------------------------------------------------------------------------
    // Map DOFs back to sample points using DOF coordinates
    // -------------------------------------------------------------------------
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
    
    std::cout << "[InterpolationManager] Extracted " << values.size() 
              << " samples along line '" << profile.name 
              << "' (length = " << total_length 
              << ", z-range = [" << M_z_min << ", " << M_z_max << "])" << std::endl;
}

// ============================================================================
// GMM Helper: Solve Normal Equations (replaces Eigen QR)
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
// Approach 1: Fit Polynomial Using OLS (GMM version)
// ============================================================================
PolynomialFit InterpolationManager::fitPolynomial(const std::vector<scalar_type>& arc_coords,
                                                  const std::vector<scalar_type>& values,
                                                  size_type order) {
    PolynomialFit result;
    size_type n = values.size();
    size_type m = order + 1;  // Number of coefficients
    
    if (n < m) {
        std::cerr << "[InterpolationManager] Error: Not enough points (" << n 
                  << ") for polynomial order " << order << std::endl;
        return result;
    }
    
    // Normalize arc coordinates to [0, 1]
    scalar_type max_arc = arc_coords.back();
    result.arc_length = max_arc;
    if (max_arc < 1.0e-15) max_arc = 1.0;  // Avoid division by zero
    
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
    // y_pred = V * coeffs
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
    
    std::cout << "[InterpolationManager] Polynomial fit (order " << order 
              << "): R² = " << result.r_squared 
              << ", coefficients = [";
    for (size_type i = 0; i < result.coefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << result.coefficients[i];
    }
    std::cout << "]" << std::endl;
    
    return result;
}

// ============================================================================
// Approach 1: Apply Polynomial as BC on Target Mesh
// ============================================================================
void InterpolationManager::applyPolynomialBC(const PolynomialFit& fit,
                                             const LineProfile& profile,
                                             const getfem::mesh_fem& mf_target,
                                             scalarVectorPtr_Type& bc_values) {
    size_type nb_dof = mf_target.nb_dof();
    bc_values.reset(new scalarVector_Type(nb_dof, 0.0));  // Initialize to zero
    
    // For cylindrical geometry, we assume the line is vertical (along z-axis)
    // So we just use the z-coordinate for projection
    
    // Find min and max z from the line profile
    scalar_type line_z_min = std::min(profile.start_point[2], profile.end_point[2]);
    scalar_type line_z_max = std::max(profile.start_point[2], profile.end_point[2]);
    
    // For each DOF, evaluate polynomial based on its z-coordinate
    for (size_type i = 0; i < nb_dof; ++i) {
        bgeot::base_node dof_pt = mf_target.point_of_basic_dof(i);
        scalar_type z = dof_pt[2];
        
        // Normalize z to [0, 1] range of the line
        scalar_type t = 0.0;
        if (std::abs(line_z_max - line_z_min) > 1e-15) {
            t = (z - line_z_min) / (line_z_max - line_z_min);
        }
        t = std::max(0.0, std::min(1.0, t));  // Clamp to [0, 1]
        
        // Evaluate polynomial
        (*bc_values)[i] = fit.evaluate(t);
    }
}

// ============================================================================
// Approach 1: Full Line Interpolation Pipeline
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
    
    // For now, use the first line profile
    const LineProfile& profile = M_lineProfiles[0];
    
    std::cout << "[InterpolationManager] Starting line interpolation via '" 
              << profile.name << "'" << std::endl;
    
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
    
    std::cout << "[InterpolationManager] Line interpolation complete. "
              << "Target solution has " << target_solution->size() << " DOFs" << std::endl;
}

// ============================================================================
// Approach 2: Direct Mesh-to-Mesh Interpolation
// ============================================================================
void InterpolationManager::interpolateViaMesh(const scalarVectorPtr_Type& source_solution,
                                              const getfem::mesh_fem& mf_source,
                                              scalarVectorPtr_Type& target_solution,
                                              const getfem::mesh_fem& mf_target) {
    // Allocate target vector
    size_type nb_dof_target = mf_target.nb_dof();
    target_solution.reset(new scalarVector_Type(nb_dof_target, 0.0));
    
    // Use GetFEM's built-in interpolation
    try {
        getfem::interpolation(mf_source, mf_target, *source_solution, *target_solution);
        std::cout << "[InterpolationManager] Mesh interpolation: " 
                  << mf_source.nb_dof() << " DOFs -> " 
                  << nb_dof_target << " DOFs" << std::endl;
    } catch (std::exception& e) {
        std::cerr << "[InterpolationManager] Mesh interpolation failed: " << e.what() << std::endl;
    }
}

// ============================================================================
// Approach 2: Build Interpolation Matrix
// ============================================================================
void InterpolationManager::buildInterpolationMatrix(const getfem::mesh_fem& mf_source,
                                                    const getfem::mesh_fem& mf_target) {
    size_type nb_dof_source = mf_source.nb_dof();
    size_type nb_dof_target = mf_target.nb_dof();
    
    // Resize interpolation matrix
    gmm::resize(M_interpMatrix, nb_dof_target, nb_dof_source);
    gmm::clear(M_interpMatrix);
    
    // Build the interpolation matrix using GetFEM
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

// ============================================================================
// Approach 2: Fast Interpolation Using Pre-computed Matrix
// ============================================================================
void InterpolationManager::interpolateViaMatrix(const scalarVectorPtr_Type& source_solution,
                                                scalarVectorPtr_Type& target_solution) {
    if (!M_matrixBuilt) {
        std::cerr << "[InterpolationManager] Error: Interpolation matrix not built!" << std::endl;
        return;
    }
    
    size_type nb_dof_target = gmm::mat_nrows(M_interpMatrix);
    target_solution.reset(new scalarVector_Type(nb_dof_target));
    
    // target = M * source
    gmm::mult(M_interpMatrix, *source_solution, *target_solution);
    
    std::cout << "[InterpolationManager] Matrix interpolation: " 
              << source_solution->size() << " -> " 
              << nb_dof_target << " DOFs" << std::endl;
}

// ============================================================================
// Unified Interface
// ============================================================================
void InterpolationManager::interpolate(const scalarVectorPtr_Type& source_solution,
                                       const getfem::mesh_fem& mf_source,
                                       scalarVectorPtr_Type& target_solution,
                                       const getfem::mesh_fem& mf_target) {
    std::cout << "[InterpolationManager] Starting interpolation..." << std::endl;
    
    if (M_approach == Approach::LINE_INTERPOLATION) {
        interpolateViaLines(source_solution, mf_source, target_solution, mf_target);
    } else {
        // Use matrix if already built, otherwise build it
        if (!M_matrixBuilt) {
            buildInterpolationMatrix(mf_source, mf_target);
        }
        interpolateViaMatrix(source_solution, target_solution);
    }
}