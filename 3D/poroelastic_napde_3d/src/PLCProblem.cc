// ============================================================================
// PLCProblem.cc - Implementation of lacuno-canalicular porosity problem
// ============================================================================
// This version supports polynomial-based boundary conditions from RussianDoll
// coupling. The outer wall BC can be set via polynomial coefficients p_v(z).
// ============================================================================

#include "../include/PLCProblem.h"
#include <iostream>
#include <cmath>

// ============================================================================
// Constructor
// ============================================================================
PLCProblem::PLCProblem(const GetPot& dataFile, 
                       Bulk* bulk, 
                       TimeLoop* time,
                       const std::string& section)
    : CoupledProblem(dataFile, bulk, time, section),  // Call base class constructor
      M_gamma(0.0),
      M_couplingSource(nullptr),
      M_pressureMass(nullptr),
      M_pressureMassBuilt(false),
      M_outerWallBCCallback(nullptr),
      M_z_min(0.0),
      M_z_max(1.0),
      M_outerWallRegion(0),  // Typically region 0 is outer wall
      M_usePolynomialBC(false)
{
    // Read leakage coefficient from data file
    M_gamma = dataFile((section + "leakage_coefficient").c_str(), 
                       dataFile("leakage_coefficient", 0.0));
    
    // Read outer wall region ID
    M_outerWallRegion = dataFile((section + "outer_wall_region").c_str(), 0);
    
    std::cout << "=== PLCProblem Created ===" << std::endl;
    std::cout << "  Section: " << section << std::endl;
    std::cout << "  Leakage coefficient (gamma): " << M_gamma << std::endl;
    std::cout << "  Outer wall region: " << M_outerWallRegion << std::endl;
}

// ============================================================================
// Boundary Condition Methods
// ============================================================================
void PLCProblem::setOuterWallBCCallback(BCCallback callback) {
    M_outerWallBCCallback = callback;
    M_usePolynomialBC = (callback != nullptr);
    
    std::cout << "[PLCProblem] Outer wall BC callback " 
              << (M_usePolynomialBC ? "set" : "cleared") << std::endl;
}

void PLCProblem::setOuterWallBCCoefficients(const std::vector<scalar_type>& coefficients,
                                            scalar_type z_min, scalar_type z_max) {
    M_bcCoefficients = coefficients;
    M_z_min = z_min;
    M_z_max = z_max;
    
    // Create internal callback that uses these coefficients
    M_outerWallBCCallback = [this](scalar_type z) -> scalar_type {
        return this->evaluateOuterWallBC(z);
    };
    
    M_usePolynomialBC = true;
    
    std::cout << "[PLCProblem] Outer wall BC coefficients set:" << std::endl;
    std::cout << "  Z range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    std::cout << "  Coefficients: [";
    for (size_t i = 0; i < M_bcCoefficients.size(); ++i) {
        std::cout << M_bcCoefficients[i];
        if (i < M_bcCoefficients.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

scalar_type PLCProblem::evaluateOuterWallBC(scalar_type z) const {
    if (M_bcCoefficients.empty()) {
        return 0.0;
    }
    
    // Normalize z to [0, 1]
    scalar_type t = 0.0;
    if (std::abs(M_z_max - M_z_min) > 1e-15) {
        t = (z - M_z_min) / (M_z_max - M_z_min);
    }
    
    // Clamp to [0, 1]
    t = std::max(0.0, std::min(1.0, t));
    
    // Evaluate polynomial: p(t) = c0 + c1*t + c2*t^2 + ...
    scalar_type result = 0.0;
    scalar_type t_power = 1.0;
    
    for (size_type i = 0; i < M_bcCoefficients.size(); ++i) {
        result += M_bcCoefficients[i] * t_power;
        t_power *= t;
    }
    
    return result;
}

// ============================================================================
// Coupling Methods (for optional two-way coupling)
// ============================================================================
void PLCProblem::setLeakageCoefficient(scalar_type gamma) {
    M_gamma = gamma;
    std::cout << "[PLCProblem] Leakage coefficient set to: " << M_gamma << std::endl;
}

void PLCProblem::setCouplingSource(const scalarVectorPtr_Type& pv_on_PLC) {
    if (!pv_on_PLC) {
        std::cerr << "[PLCProblem] Warning: Null coupling source provided!" << std::endl;
        return;
    }
    
    size_type expectedSize = getNbPressureDOF();
    if (pv_on_PLC->size() != expectedSize) {
        std::cerr << "[PLCProblem] Warning: Coupling source size mismatch! "
                  << "Expected " << expectedSize 
                  << ", got " << pv_on_PLC->size() << std::endl;
    }
    
    // Copy the coupling source
    M_couplingSource.reset(new scalarVector_Type(*pv_on_PLC));    
    std::cout << "[PLCProblem] Coupling source set (norm = " 
              << gmm::vect_norm2(*M_couplingSource) << ")" << std::endl;
}

// ============================================================================
// Build Pressure Mass Matrix
// ============================================================================
void PLCProblem::buildPressureMassMatrix() {
    if (M_pressureMassBuilt) {
        return;  // Already built
    }
    
    if (!M_DarcyPB) {
        std::cerr << "[PLCProblem] Error: Darcy problem not set!" << std::endl;
        return;
    }
    
    size_type nbPressureDOF = getNbPressureDOF();
    
    // Allocate pressure mass matrix
    M_pressureMass.reset(new sparseMatrix_Type(nbPressureDOF, nbPressureDOF));
    gmm::clear(*M_pressureMass);
    
    // Get pressure FEM and integration method
    FEM* pressureFEM = M_DarcyPB->getFEM("Pressure");
    
    // Build mass matrix: M_ij = ∫ N_i * N_j dV
    massL2Standard(M_pressureMass, 
                   *pressureFEM,    // Test FEM
                   *pressureFEM,    // Trial FEM (same for mass matrix)
                   M_intMethod);
    
    M_pressureMassBuilt = true;
    
    std::cout << "[PLCProblem] Pressure mass matrix built (" 
              << nbPressureDOF << " x " << nbPressureDOF << ")" << std::endl;
}

// ============================================================================
// Assemble Matrix (Override)
// ============================================================================
void PLCProblem::assembleMatrix() {
    // Call base class to assemble standard poroelasticity matrices
    CoupledProblem::assembleMatrix();
    
    // Also build the pressure mass matrix for potential coupling term
    if (std::abs(M_gamma) > 1e-15) {
        buildPressureMassMatrix();
    }
}

// ============================================================================
// Assemble RHS (Override)
// ============================================================================
void PLCProblem::assembleRHS() {
    // Call base class to assemble standard RHS
    CoupledProblem::assembleRHS();
    
    // Add the inter-porosity coupling term (for two-way coupling, if enabled)
    if (std::abs(M_gamma) > 1e-15 && hasCouplingSource()) {
        assembleCouplingRHS();
    }
}

// ============================================================================
// Assemble Coupling RHS Term
// ============================================================================
void PLCProblem::assembleCouplingRHS() {
    // Add coupling term: +γ * M * p_v to the pressure equation RHS
    // This is for two-way coupling (typically not used in one-way mode)
    
    if (!hasCouplingSource()) {
        return;
    }
    
    if (std::abs(M_gamma) < 1e-15) {
        return;  // Zero leakage coefficient - skip
    }
    
    if (!M_pressureMassBuilt) {
        std::cerr << "[PLCProblem] Warning: Pressure mass matrix not built! "
                  << "Call assembleMatrix() first." << std::endl;
        buildPressureMassMatrix();
    }
    
    size_type nbPressureDOF = getNbPressureDOF();
    size_type nbVelocityDOF = getNbVelocityDOF();
    size_type nbElastDOF = getNbElastDOF();
    
    // Compute coupling RHS: γ * M * p_v
    scalarVectorPtr_Type couplingRHS;
    couplingRHS.reset(new scalarVector_Type(nbPressureDOF, 0.0));
    
    // Matrix-vector multiplication: couplingRHS = M * p_v
    gmm::mult(*M_pressureMass, *M_couplingSource, *couplingRHS);
    
    // Scale by leakage coefficient: couplingRHS = γ * M * p_v
    gmm::scale(*couplingRHS, M_gamma);
    
    // Add to global RHS at the correct position
    // DOFs are ordered as: [Elasticity | Velocity | Pressure]
    size_type pressureOffset = nbElastDOF + nbVelocityDOF;
    
    scalarVectorPtr_Type globalRHS = M_Sys->getRHS();
    for (size_type i = 0; i < nbPressureDOF; ++i) {
        (*globalRHS)[pressureOffset + i] += (*couplingRHS)[i];
    }
    
    std::cout << "[PLCProblem] Coupling RHS assembled (norm = " 
              << gmm::vect_norm2(*couplingRHS) << ")" << std::endl;
}

// ============================================================================
// Enforce Boundary Conditions (Override)
// ============================================================================
void PLCProblem::enforceStrongBC(bool firstTime) {
    // First, enforce standard BCs from base class (elasticity)
    if (firstTime) {
        M_ElastPB->enforceStrongBC(true);
    } else {
        M_ElastPB->enforceStrongBC(false);
    }
    
    // Then enforce polynomial BC on outer wall for Darcy pressure
    if (M_usePolynomialBC) {
        enforceOuterWallBC(firstTime);
    }
}

// ============================================================================
// Enforce Polynomial Dirichlet BC on Outer Wall
// ============================================================================
void PLCProblem::enforceOuterWallBC(bool firstTime) {
    std::cout << "[PLCProblem] Enforcing polynomial BC on outer wall (region " 
              << M_outerWallRegion << ")..." << std::endl;
    
    if (!M_outerWallBCCallback) {
        std::cerr << "[PLCProblem] Warning: No BC callback set for outer wall!" << std::endl;
        return;
    }
    
    if (!M_DarcyPB) {
        std::cerr << "[PLCProblem] Error: Darcy problem not set!" << std::endl;
        return;
    }
    
    // Get pressure FEM
    const getfem::mesh_fem& mf_pressure = *(M_DarcyPB->getFEM("Pressure")->getFEM());
    
    // Get system matrix and RHS
    sparseMatrixPtr_Type matrix = M_Sys->getMatrix();
    scalarVectorPtr_Type rhs = M_Sys->getRHS();
    
    // Calculate offset for pressure DOFs in global system
    // DOFs are ordered as: [Elasticity | Velocity | Pressure]
    size_type nbElastDOF = getNbElastDOF();
    size_type nbVelocityDOF = getNbVelocityDOF();
    size_type pressureOffset = nbElastDOF + nbVelocityDOF;
    
    // Get mesh and boundary region
    const getfem::mesh& mesh = mf_pressure.linked_mesh();
    getfem::mesh_region region = mesh.region(M_outerWallRegion);
    
    // Find all pressure DOFs on the outer wall boundary
    dal::bit_vector boundary_dofs;
    
    // Iterate over faces in the outer wall region
    for (getfem::mr_visitor it(region); !it.finished(); ++it) {
        if (!it.is_face()) continue;
        
        size_type cv = it.cv();
        short int f = it.f();
        
        // Get DOFs on this face
        getfem::pfem pf = mf_pressure.fem_of_element(cv);
        if (pf == nullptr) continue;
        
        // Get local DOF indices on this face
        bgeot::convex<base_node> cv_ref = pf->node_convex(cv);
        
        // Get DOF indices for this element
        auto dof_indices = mf_pressure.ind_basic_dof_of_element(cv);
        
        // Iterate over DOFs and check if they're on the boundary face
        for (size_type i = 0; i < dof_indices.size(); ++i) {
            size_type dof = dof_indices[i];
            bgeot::base_node dof_pt = mf_pressure.point_of_basic_dof(dof);
            
            // Check if this DOF is on the face (approximately)
            // For pressure DOFs on the boundary, we include all DOFs whose
            // points are on the outer surface
            
            // Simple check: if the DOF point is on the outer wall
            // This depends on geometry - for a cylinder, check radius
            // For generality, we mark all DOFs on faces in the region
            
            boundary_dofs.add(dof);
        }
    }
    
    // Apply Dirichlet BC to each boundary DOF
    size_type bc_count = 0;
    
    for (dal::bv_visitor dof(boundary_dofs); !dof.finished(); ++dof) {
        size_type global_dof = pressureOffset + dof;
        
        // Get DOF coordinates
        bgeot::base_node pt = mf_pressure.point_of_basic_dof(dof);
        scalar_type z = pt[2];  // Assuming z is the axial direction
        
        // Evaluate polynomial BC
        scalar_type bc_value = M_outerWallBCCallback(z);
        
        // Enforce Dirichlet BC using penalty or elimination method
        if (firstTime) {
            // Modify matrix: set row to identity
            size_type ncols = gmm::mat_ncols(*matrix);
            
            // Clear row (except diagonal)
            for (size_type j = 0; j < ncols; ++j) {
                (*matrix)(global_dof, j) = 0.0;
            }
            
            // Set diagonal to 1
            (*matrix)(global_dof, global_dof) = 1.0;
        }
        
        // Set RHS
        (*rhs)[global_dof] = bc_value;
        
        bc_count++;
    }
    
    std::cout << "[PLCProblem] Applied polynomial BC to " << bc_count 
              << " DOFs on outer wall" << std::endl;
    
    // Sample BC values for debugging
    if (bc_count > 0) {
        scalar_type p_min = M_outerWallBCCallback(M_z_min);
        scalar_type p_max = M_outerWallBCCallback(M_z_max);
        scalar_type p_mid = M_outerWallBCCallback((M_z_min + M_z_max) / 2.0);
        
        std::cout << "  BC values: p(z_min)=" << p_min 
                  << ", p(z_mid)=" << p_mid 
                  << ", p(z_max)=" << p_max << std::endl;
    }
}