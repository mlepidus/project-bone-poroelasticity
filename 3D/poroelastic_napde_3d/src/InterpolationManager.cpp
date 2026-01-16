// ============================================================================
// InterpolationManager.cpp - Implementation of solution transfer between meshes
// ============================================================================

#include "InterpolationManager.h"
#include <getfem/getfem_interpolation.h>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <iostream>
#include <cmath>

// ============================================================================
// Constructor
// ============================================================================
InterpolationManager::InterpolationManager(const GetPot& dataFile,
                                           Bulk* sourceBulk,
                                           Bulk* targetBulk)
    : M_sourceBulk(sourceBulk),
      M_targetBulk(targetBulk),
      M_matrixBuilt(false),
      M_section("russian_doll/interpolation/")
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

void InterpolationManager::addLineProfile(const LineProfile& profile) {
    M_lineProfiles.push_back(profile);
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
// Approach 1: Extract Solution Along a Line
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
    
    // Compute total arc length
    bgeot::base_node diff = profile.end_point - profile.start_point;
    scalar_type total_length = gmm::vect_norm2(diff);
    
    // Create temporary 1D mesh for interpolation target
    getfem::mesh mesh_1d;
    std::vector<size_type> point_indices;
    
    for (size_type i = 0; i < profile.num_samples; ++i) {
        scalar_type t = static_cast<scalar_type>(i) / (profile.num_samples - 1);
        bgeot::base_node pt = interpolatePoint(profile, t);
        point_indices.push_back(mesh_1d.add_point(pt));
        arc_coords.push_back(t * total_length);
    }
    
    // Add 1D elements (optional, for mesh_fem)
    for (size_type i = 0; i < profile.num_samples - 1; ++i) {
        std::vector<size_type> ind = {point_indices[i], point_indices[i+1]};
        mesh_1d.add_convex(bgeot::simplex_geotrans(1, 1), ind.begin());
    }
    
    // Create mesh_fem on 1D mesh
    getfem::mesh_fem mf_1d(mesh_1d, 1);
    mf_1d.set_classical_finite_element(1);  // P1 elements
    
    // Allocate target vector
    std::vector<scalar_type> interp_values(mf_1d.nb_dof());
    
    // Perform interpolation from 3D to 1D
    getfem::interpolation(mf_source, mf_1d, *solution, interp_values);
    
    // Extract values at nodes (P1 DOFs are at vertices)
    values.resize(profile.num_samples);
    for (size_type i = 0; i < profile.num_samples; ++i) {
        // For P1 elements, DOF i corresponds to vertex i
        values[i] = interp_values[i];
    }
    
    std::cout << "[InterpolationManager] Extracted " << values.size() 
              << " samples along line '" << profile.name << "'" << std::endl;
}

// ============================================================================
// Approach 1: Fit Polynomial Using OLS
// ============================================================================
PolynomialFit InterpolationManager::fitPolynomial(const std::vector<scalar_type>& arc_coords,
                                                  const std::vector<scalar_type>& values,
                                                  size_type order) {
    PolynomialFit result;
    size_type n = values.size();
    
    // Normalize arc coordinates to [0, 1]
    scalar_type max_arc = arc_coords.back();
    result.arc_length = max_arc;
    
    // Build Vandermonde matrix for polynomial fitting
    Eigen::MatrixXd V(n, order + 1);
    Eigen::VectorXd y(n);
    
    for (size_type i = 0; i < n; ++i) {
        scalar_type t = arc_coords[i] / max_arc;  // Normalized coordinate
        scalar_type t_power = 1.0;
        for (size_type j = 0; j <= order; ++j) {
            V(i, j) = t_power;
            t_power *= t;
        }
        y(i) = values[i];
    }
    
    // Solve least squares: V * coeffs = y
    // Using QR decomposition for stability
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(V);
    result.coefficients = qr.solve(y);
    
    // Compute R² goodness of fit
    Eigen::VectorXd y_pred = V * result.coefficients;
    Eigen::VectorXd residuals = y - y_pred;
    
    scalar_type ss_res = residuals.squaredNorm();
    scalar_type mean_y = y.mean();
    scalar_type ss_tot = (y.array() - mean_y).matrix().squaredNorm();
    
    result.r_squared = 1.0 - ss_res / ss_tot;
    
    std::cout << "[InterpolationManager] Polynomial fit (order " << order 
              << "): R² = " << result.r_squared << std::endl;
    
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
    bc_values = std::make_shared<scalarVector_Type>(nb_dof);
    
    // Compute direction vector of the line
    bgeot::base_node dir = profile.end_point - profile.start_point;
    scalar_type line_length = gmm::vect_norm2(dir);
    
    // Normalize direction
    for (int d = 0; d < 3; ++d) {
        dir[d] /= line_length;
    }
    
    // For each DOF, find its projection onto the line and evaluate polynomial
    for (size_type i = 0; i < nb_dof; ++i) {
        bgeot::base_node dof_pt = mf_target.point_of_basic_dof(i);
        
        // Project DOF point onto line
        bgeot::base_node v = dof_pt - profile.start_point;
        scalar_type projection = gmm::vect_sp(v, dir);  // Dot product
        
        // Clamp to [0, line_length] and normalize
        scalar_type t = std::max(0.0, std::min(projection, line_length)) / line_length;
        
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
    // TODO: Support multiple lines with averaging or selection logic
    const LineProfile& profile = M_lineProfiles[0];
    
    // Step 1: Extract values along line
    std::vector<scalar_type> arc_coords, values;
    extractAlongLine(source_solution, mf_source, profile, arc_coords, values);
    
    // Step 2: Fit polynomial
    PolynomialFit fit = fitPolynomial(arc_coords, values, profile.polynomial_order);
    M_polynomialFits.push_back(fit);
    
    // Step 3: Apply to target mesh
    applyPolynomialBC(fit, profile, mf_target, target_solution);
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
    target_solution = std::make_shared<scalarVector_Type>(nb_dof_target);
    
    // Use GetFEM's built-in interpolation
    // This handles the case where target mesh is inside source mesh
    getfem::interpolation(mf_source, mf_target, *source_solution, *target_solution);
    
    std::cout << "[InterpolationManager] Mesh interpolation: " 
              << mf_source.nb_dof() << " DOFs -> " 
              << nb_dof_target << " DOFs" << std::endl;
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
    // This computes M such that: target = M * source
    getfem::interpolation(mf_source, mf_target, M_interpMatrix);
    
    M_matrixBuilt = true;
    
    std::cout << "[InterpolationManager] Built interpolation matrix: " 
              << nb_dof_target << " x " << nb_dof_source << std::endl;
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
    target_solution = std::make_shared<scalarVector_Type>(nb_dof_target);
    
    // target = M * source
    gmm::mult(M_interpMatrix, *source_solution, *target_solution);
}

// ============================================================================
// Unified Interface
// ============================================================================
void InterpolationManager::interpolate(const scalarVectorPtr_Type& source_solution,
                                       const getfem::mesh_fem& mf_source,
                                       scalarVectorPtr_Type& target_solution,
                                       const getfem::mesh_fem& mf_target) {
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