// ============================================================================
// PLCProblem.cc - Implementation of lacuno-canalicular porosity problem
// ============================================================================
// Enhanced version with:
// - Polynomial BC from PV pressure (outer wall Dirichlet)
// - Coupling RHS term gamma * M * p_v for mass balance equation
// - Coupling matrix term +gamma * M added to system
// ============================================================================

#include "../include/PLCProblem.h"
#include <iostream>
#include <cmath>
#include <algorithm>

// ============================================================================
// Constructor
// ============================================================================
PLCProblem::PLCProblem(const GetPot& dataFile, 
                       Bulk* bulk, 
                       TimeLoop* time,
                       const std::string& section)
    : CoupledProblem(dataFile, bulk, time, section),
      M_gamma(0.0),
      M_couplingSource(nullptr),
      M_pressureMass(nullptr),
      M_pressureMassBuilt(false),
      M_couplingMatrixAdded(false),
      M_outerWallBCCallbackP(nullptr),
      M_outerWallBCCallbackX(nullptr),
      M_outerWallBCCallbackY(nullptr),
      M_outerWallBCCallbackZ(nullptr),
      M_z_min(0.0),
      M_z_max(1.0),
      M_outerWallRegion(0),
      M_usePolynomialBC(false)
{
    // Read leakage coefficient from data file
    // Try section-specific first, then global
    M_gamma = dataFile((section + "darcy/leakage").c_str(), 
                       dataFile("leakage_coefficient", 0.0));
    
    // Read outer wall region ID
    M_outerWallRegion = dataFile((section + "outer_wall_region").c_str(), 0);
    
    std::cout << "=====================================================" << std::endl;
    std::cout << " PLCProblem Created" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Section: " << section << std::endl;
    std::cout << "  Leakage coefficient (gamma): " << M_gamma << std::endl;
    std::cout << "  Outer wall region: " << M_outerWallRegion << std::endl;
    std::cout << "=====================================================" << std::endl;
}

// ============================================================================
// Boundary Condition Methods
// ============================================================================
void PLCProblem::setOuterWallBCCallback(BCCallback callback) {
    M_outerWallBCCallbackP = callback;
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
    M_outerWallBCCallbackP = [this](scalar_type z) -> scalar_type {
        return this->evaluateOuterWallBC(z);
    };
    
    M_usePolynomialBC = true;
    
    std::cout << "[PLCProblem] Outer wall BC coefficients set:" << std::endl;
    std::cout << "  Z range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    std::cout << "  Polynomial order: " << (M_bcCoefficients.size() - 1) << std::endl;
    std::cout << "  Coefficients: [";
    for (size_t i = 0; i < M_bcCoefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << M_bcCoefficients[i];
    }
    std::cout << "]" << std::endl;
    
    // Also update the coupling source using these coefficients
    setCouplingSourceFromPolynomial(coefficients, z_min, z_max);
}

// ============================================================================
// Set Displacement Boundary Conditions
// ============================================================================
void PLCProblem::setOuterWallDisplacementBCCoefficients(
    const std::vector<scalar_type>& coeffs_x,
    const std::vector<scalar_type>& coeffs_y,
    const std::vector<scalar_type>& coeffs_z,
    scalar_type z_min, scalar_type z_max) {
    
    M_dispXCoefficients = coeffs_x;
    M_dispYCoefficients = coeffs_y;
    M_dispZCoefficients = coeffs_z;
    M_z_min = z_min;  // Assuming same z-range for displacement
    M_z_max = z_max;
    
    // Create internal callbacks that use these coefficients
    M_outerWallBCCallbackX = [this](scalar_type z) -> scalar_type {
        return this->evaluatePolynomial(z, M_dispXCoefficients, M_z_min, M_z_max);
    };
    M_outerWallBCCallbackY = [this](scalar_type z) -> scalar_type {
        return this->evaluatePolynomial(z, M_dispYCoefficients, M_z_min, M_z_max);
    };
    M_outerWallBCCallbackZ = [this](scalar_type z) -> scalar_type {
        return this->evaluatePolynomial(z, M_dispZCoefficients, M_z_min, M_z_max);
    };
    
    M_useDisplacementBC = true;
    
    std::cout << "[PLCProblem] Outer wall displacement BC coefficients set:" << std::endl;
    std::cout << "  Z range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    std::cout << "  Polynomial order: " << (M_dispXCoefficients.size() - 1) << std::endl;
    
    // Also compute displacement at PLC DOFs if needed
    computeDisplacementOnPLCDOFs();
}

// ============================================================================
// Compute Displacement at PLC DOF Locations
// ============================================================================
void PLCProblem::computeDisplacementOnPLCDOFs() {
    // This is optional - only needed if you want to use displacement for something else
    if (!M_ElastPB) {
        std::cerr << "[PLCProblem] Error: Elasticity problem not set!" << std::endl;
        return;
    }
    
    // Get displacement mesh_fem
    const getfem::mesh_fem& mf_displacement = *(M_ElastPB->getFEM()->getFEM());
    size_type nbDisplacementDOF = mf_displacement.nb_dof();
    
    // Note: This would store displacement values if needed for coupling
    // Currently not used but available if needed
}

// ============================================================================
// Enforce Polynomial Displacement BC on Outer Wall
// ============================================================================
void PLCProblem::enforceOuterWallDisplacementBC(bool firstTime) {
    std::cout << "[PLCProblem] Enforcing polynomial displacement BC on outer wall (region " 
              << M_outerWallRegion << ")..." << std::endl;
    
    if (!M_outerWallBCCallbackX && !M_outerWallBCCallbackY && !M_outerWallBCCallbackZ) {
        std::cerr << "[PLCProblem] Warning: No displacement BC callbacks set!" << std::endl;
        return;
    }
    
    if (!M_ElastPB) {
        std::cerr << "[PLCProblem] Error: Elasticity problem not set!" << std::endl;
        return;
    }
    
    // Get displacement FEM
    const getfem::mesh_fem& mf_displacement = *(M_ElastPB->getFEM()->getFEM());
    
    // Get system matrix and RHS
    sparseMatrixPtr_Type matrix = M_Sys->getMatrix();
    scalarVectorPtr_Type rhs = M_Sys->getRHS();
    
    // Calculate offset for displacement DOFs in global system
    // DOFs are ordered as: [Elasticity | Velocity | Pressure]
    size_type dispOffset = 0;  // Displacement DOFs come first
    
    // Get mesh and boundary region
    const getfem::mesh& mesh = mf_displacement.linked_mesh();
    
    // Check if region exists
    if (!mesh.has_region(M_outerWallRegion)) {
        std::cerr << "[PLCProblem] Warning: Region " << M_outerWallRegion 
                  << " does not exist!" << std::endl;
        return;
    }
    
    getfem::mesh_region region = mesh.region(M_outerWallRegion);
    
    // Find all displacement DOFs on the outer wall boundary
    dal::bit_vector boundary_dofs;
    
    // Iterate over faces in the outer wall region
    for (getfem::mr_visitor it(region); !it.finished(); ++it) {
        if (!it.is_face()) continue;
        
        size_type cv = it.cv();
        
        // Check if element has FEM
        getfem::pfem pf = mf_displacement.fem_of_element(cv);
        if (pf == nullptr) continue;
        
        // Get DOF indices for this element
        auto dof_indices = mf_displacement.ind_basic_dof_of_element(cv);
        
        // Mark all DOFs of this element on the boundary
        for (size_type i = 0; i < dof_indices.size(); ++i) {
            boundary_dofs.add(dof_indices[i]);
        }
    }
    
    // Apply Dirichlet BC to each boundary DOF
    size_type bc_count = 0;
    size_type ncols = gmm::mat_ncols(*matrix);
    
    for (dal::bv_visitor dof(boundary_dofs); !dof.finished(); ++dof) {
        // Get DOF coordinates
        bgeot::base_node pt = mf_displacement.point_of_basic_dof(dof);
        scalar_type z = pt[2];  // Assuming z is the axial direction
        
        // Get the component of this DOF (0=x, 1=y, 2=z)
        // This depends on how your FEM is set up
        size_type component = dof % 3;  // Assuming interleaved: ux1, uy1, uz1, ux2, uy2, uz2, ...
        size_type global_dof = dispOffset + dof;
        
        // Evaluate polynomial BC for the appropriate component
        scalar_type bc_value = 0.0;
        switch (component) {
            case 0:  // x-component
                if (M_outerWallBCCallbackX) {
                    bc_value = M_outerWallBCCallbackX(z);
                }
                break;
            case 1:  // y-component
                if (M_outerWallBCCallbackY) {
                    bc_value = M_outerWallBCCallbackY(z);
                }
                break;
            case 2:  // z-component
                if (M_outerWallBCCallbackZ) {
                    bc_value = M_outerWallBCCallbackZ(z);
                }
                break;
        }
        
        // Enforce Dirichlet BC
        if (firstTime) {
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
    
    std::cout << "[PLCProblem] Applied displacement polynomial BC to " << bc_count 
              << " DOFs on outer wall" << std::endl;
    
    // Sample BC values for debugging
    if (bc_count > 0) {
        scalar_type ux_min = M_outerWallBCCallbackX ? M_outerWallBCCallbackX(M_z_min) : 0.0;
        scalar_type uy_min = M_outerWallBCCallbackY ? M_outerWallBCCallbackY(M_z_min) : 0.0;
        scalar_type uz_min = M_outerWallBCCallbackZ ? M_outerWallBCCallbackZ(M_z_min) : 0.0;
        
        scalar_type ux_mid = M_outerWallBCCallbackX ? 
            M_outerWallBCCallbackX((M_z_min + M_z_max) / 2.0) : 0.0;
        scalar_type uy_mid = M_outerWallBCCallbackY ? 
            M_outerWallBCCallbackY((M_z_min + M_z_max) / 2.0) : 0.0;
        scalar_type uz_mid = M_outerWallBCCallbackZ ? 
            M_outerWallBCCallbackZ((M_z_min + M_z_max) / 2.0) : 0.0;
        
        scalar_type ux_max = M_outerWallBCCallbackX ? M_outerWallBCCallbackX(M_z_max) : 0.0;
        scalar_type uy_max = M_outerWallBCCallbackY ? M_outerWallBCCallbackY(M_z_max) : 0.0;
        scalar_type uz_max = M_outerWallBCCallbackZ ? M_outerWallBCCallbackZ(M_z_max) : 0.0;
        
        std::cout << "  Displacement BC sample values:" << std::endl;
        std::cout << "    u(z_min=" << M_z_min << "): [" 
                  << ux_min << ", " << uy_min << ", " << uz_min << "]" << std::endl;
        std::cout << "    u(z_mid=" << (M_z_min + M_z_max)/2.0 << "): [" 
                  << ux_mid << ", " << uy_mid << ", " << uz_mid << "]" << std::endl;
        std::cout << "    u(z_max=" << M_z_max << "): [" 
                  << ux_max << ", " << uy_max << ", " << uz_max << "]" << std::endl;
    }
}

scalar_type PLCProblem::evaluateOuterWallBC(scalar_type z) const {
    return evaluatePolynomial(z, M_bcCoefficients, M_z_min, M_z_max);
}

scalar_type PLCProblem::evaluatePolynomial(scalar_type z,
                                           const std::vector<scalar_type>& coeffs,
                                           scalar_type z_min,
                                           scalar_type z_max) const {
    if (coeffs.empty()) {
        return 0.0;
    }
    
    // Normalize z to [0, 1]
    scalar_type t = 0.0;
    if (std::abs(z_max - z_min) > 1e-15) {
        t = (z - z_min) / (z_max - z_min);
    }
    
    // Clamp to [0, 1]
    t = std::max(0.0, std::min(1.0, t));
    
    // Evaluate polynomial: p(t) = c0 + c1*t + c2*t^2 + ...
    scalar_type result = 0.0;
    scalar_type t_power = 1.0;
    
    for (size_type i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * t_power;
        t_power *= t;
    }
    
    return result;
}

// ============================================================================
// Coupling Methods
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

void PLCProblem::setCouplingSourceFromPolynomial(const std::vector<scalar_type>& coefficients,
                                                  scalar_type z_min, scalar_type z_max) {
    if (!M_DarcyPB) {
        std::cerr << "[PLCProblem] Error: Darcy problem not set!" << std::endl;
        return;
    }
    
    size_type nbPressureDOF = getNbPressureDOF();
    
    // Allocate coupling source vector
    M_couplingSource.reset(new scalarVector_Type(nbPressureDOF, 0.0));
    
    // Get pressure FEM
    const getfem::mesh_fem& mf_pressure = *(M_DarcyPB->getFEM("Pressure")->getFEM());
    
    // Evaluate polynomial at each pressure DOF location
    for (size_type i = 0; i < nbPressureDOF; ++i) {
        bgeot::base_node pt = mf_pressure.point_of_basic_dof(i);
        scalar_type z = pt[2];  // Assuming z is axial direction
        (*M_couplingSource)[i] = evaluatePolynomial(z, coefficients, z_min, z_max);
    }
    
    std::cout << "[PLCProblem] Coupling source computed from polynomial (norm = " 
              << gmm::vect_norm2(*M_couplingSource) << ")" << std::endl;
}

// ============================================================================
// Build Pressure Mass Matrix
// ============================================================================
void PLCProblem::buildPressureMassMatrix() {
    if (M_pressureMassBuilt) {
        return;
    }
    
    if (!M_DarcyPB) {
        std::cerr << "[PLCProblem] Error: Darcy problem not set!" << std::endl;
        return;
    }
    
    size_type nbPressureDOF = getNbPressureDOF();
    
    // Allocate pressure mass matrix
    M_pressureMass.reset(new sparseMatrix_Type(nbPressureDOF, nbPressureDOF));
    gmm::clear(*M_pressureMass);
    
    // Get pressure FEM
    FEM* pressureFEM = M_DarcyPB->getFEM("Pressure");
    
    // Build mass matrix: M_ij = integral N_i * N_j dV
    // Use the massL2Standard function from UsefulFunctions
    massL2Standard(M_pressureMass, 
                   *pressureFEM,
                   *pressureFEM,
                   M_intMethod);
    
    M_pressureMassBuilt = true;
    
    std::cout << "[PLCProblem] Pressure mass matrix built (" 
              << nbPressureDOF << " x " << nbPressureDOF << ")" << std::endl;
}

// ============================================================================
// Add Coupling Matrix
// ============================================================================
void PLCProblem::addCouplingMatrix() {
    if (std::abs(M_gamma) < 1e-15) {
        return;  // No coupling
    }
    
    if (!M_pressureMassBuilt) {
        buildPressureMassMatrix();
    }
    
    if (M_couplingMatrixAdded) {
        return;  // Already added
    }
    
    // Get system matrix
    sparseMatrixPtr_Type matrix = M_Sys->getMatrix();
    
    // Get DOF offsets
    size_type nbElastDOF = getNbElastDOF();
    size_type nbVelocityDOF = getNbVelocityDOF();
    size_type nbPressureDOF = getNbPressureDOF();
    size_type pressureOffset = nbElastDOF + nbVelocityDOF;
    
    // Add +gamma * M to the (pressure, pressure) block
    // This corresponds to the +gamma * M * p_l term in the mass equation
    for (size_type i = 0; i < nbPressureDOF; ++i) {
        for (size_type j = 0; j < nbPressureDOF; ++j) {
            scalar_type M_ij = (*M_pressureMass)(i, j);
            if (std::abs(M_ij) > 1e-15) {
                (*matrix)(pressureOffset + i, pressureOffset + j) += M_gamma * M_ij;
            }
        }
    }
    
    M_couplingMatrixAdded = true;
    
    std::cout << "[PLCProblem] Coupling matrix (+gamma*M) added to system" << std::endl;
}

// ============================================================================
// Assemble Matrix (Override)
// ============================================================================
void PLCProblem::assembleMatrix() {
    // Call base class to assemble standard poroelasticity matrices
    CoupledProblem::assembleMatrix();
    
    // Build the pressure mass matrix (needed for coupling)
    if (std::abs(M_gamma) > 1e-15) {
        buildPressureMassMatrix();
    }
    
    std::cout << "[PLCProblem] Matrix assembly complete" << std::endl;
}

// ============================================================================
// Assemble RHS (Override)
// ============================================================================
void PLCProblem::assembleRHS() {
    // Call base class to assemble standard RHS
    CoupledProblem::assembleRHS();
    
    // Add the inter-porosity coupling term: +gamma * M * p_v
    if (std::abs(M_gamma) > 1e-15 && hasCouplingSource()) {
        assembleCouplingRHS();
    }



//    if (M_usePolynomialBC) {
//        assembleOuterWallNeumannRHS();
//   }
}
/*
void PLCProblem::assembleOuterWallNeumannRHS() {
    if (!M_DarcyPB) {
        std::cerr << "[PLCProblem] Error: Darcy problem not set!" << std::endl;
        return;
    }
    
    std::cout << "[PLCProblem] Assembling Neumann BC on outer wall from p_v..." << std::endl;
    
    // Get velocity FEM (Neumann BC affects velocity equation)
    FEM* velocityFEM = M_DarcyPB->getFEM("Velocity");
    FEM* coeffFEM = M_DarcyPB->getFEM("Pressure");  // Using pressure FEM for coefficients
    const getfem::mesh_fem& mf_velocity = *(velocityFEM->getFEM());
    const getfem::mesh_fem& mf_coeff = *(coeffFEM->getFEM());
    
    // Get mesh and region
    const getfem::mesh& mesh = mf_velocity.linked_mesh();
    
    if (!mesh.has_region(M_outerWallRegion)) {
        std::cerr << "[PLCProblem] Warning: Outer wall region " << M_outerWallRegion 
                  << " not found!" << std::endl;
        return;
    }
    
    // Build pressure values at coefficient DOF locations
    scalarVector_Type p_bc(mf_coeff.nb_dof(), 0.0);
    for (size_type i = 0; i < mf_coeff.nb_dof(); ++i) {
        bgeot::base_node pt = mf_coeff.point_of_basic_dof(i);
        scalar_type z = pt[2];
        p_bc[i] = evaluatePolynomial(z, M_bcCoefficients, M_z_min, M_z_max);
    }
    
    // Assemble: -∫_Γ p_v * v·n dS (natural BC for pressure)
    // This is the contribution: -p_v on Neumann boundary for the mixed formulation
    
    scalarVectorPtr_Type bcVec;
    bcVec.reset(new scalarVector_Type(mf_velocity.nb_dof(), 0.0));
    
    getfem::generic_assembly assem;
    assem.set("p=data$1(#2);"
              "t=-comp(vBase(#1).Normal().Base(#2));"
              "V$1(#1)+=(t(:,i, i, k).p(k));");
    
    assem.push_mi(M_intMethod);
    assem.push_mf(mf_velocity);
    assem.push_mf(mf_coeff);
    assem.push_data(p_bc);
    assem.push_vec(*bcVec);
    
    // Assemble on outer wall region
    assem.assembly(mesh.region(M_outerWallRegion));
    
    // Add to global RHS
    // DOF order: [Elasticity | Velocity | Pressure]
    size_type nbElastDOF = getNbElastDOF();
    size_type velocityOffset = nbElastDOF;
    
    scalarVectorPtr_Type globalRHS = M_Sys->getRHS();
    for (size_type i = 0; i < mf_velocity.nb_dof(); ++i) {
        (*globalRHS)[velocityOffset + i] += (*bcVec)[i];
    }
    
    std::cout << "[PLCProblem] Outer wall Neumann BC assembled (norm = " 
              << gmm::vect_norm2(*bcVec) << ")" << std::endl;
}
*/

// ============================================================================
// Assemble Coupling RHS Term
// ============================================================================
void PLCProblem::assembleCouplingRHS() {
    // Add coupling term: +gamma * M * p_v to the pressure equation RHS
    // This represents the source term from inter-porosity exchange
    
    if (!hasCouplingSource()) {
        std::cout << "[PLCProblem] No coupling source - skipping coupling RHS" << std::endl;
        return;
    }
    
    if (std::abs(M_gamma) < 1e-15) {
        return;
    }
    
    if (!M_pressureMassBuilt) {
        std::cerr << "[PLCProblem] Warning: Pressure mass matrix not built! "
                  << "Building now..." << std::endl;
        buildPressureMassMatrix();
    }
    
    size_type nbPressureDOF = getNbPressureDOF();
    size_type nbVelocityDOF = getNbVelocityDOF();
    size_type nbElastDOF = getNbElastDOF();
    
    // Compute coupling RHS: gamma * M * p_v
    scalarVectorPtr_Type couplingRHS;
    couplingRHS.reset(new scalarVector_Type(nbPressureDOF, 0.0));
    
    // Matrix-vector multiplication: couplingRHS = M * p_v
    gmm::mult(*M_pressureMass, *M_couplingSource, *couplingRHS);
    
    // Scale by leakage coefficient: couplingRHS = gamma * M * p_v
    gmm::scale(*couplingRHS, M_gamma);
    
    // Add to global RHS at the correct position
    // DOFs are ordered as: [Elasticity | Velocity | Pressure]
    size_type pressureOffset = nbElastDOF + nbVelocityDOF;
    
    scalarVectorPtr_Type globalRHS = M_Sys->getRHS();
    for (size_type i = 0; i < nbPressureDOF; ++i) {
        (*globalRHS)[pressureOffset + i] += (*couplingRHS)[i];
    }
        // DIAGNOSTIC
    scalar_type rhs_norm = gmm::vect_norm2(*couplingRHS);
    scalar_type source_norm = gmm::vect_norm2(*M_couplingSource);
    std::cout << "[PLCProblem] Coupling RHS assembled:" << std::endl;
    std::cout << "  gamma = " << M_gamma << std::endl;
    std::cout << "  ||M * p_v|| = " << gmm::vect_norm2(*couplingRHS) / M_gamma << std::endl;
    std::cout << "  ||gamma * M * p_v|| = " << gmm::vect_norm2(*couplingRHS) << std::endl;

     // Sanity checks
    if (rhs_norm > 1e10) {
        std::cerr << "  WARNING: Coupling RHS norm is very large!" << std::endl;
    }
    if (source_norm < 1e-12) {
        std::cerr << "  WARNING: Coupling source is nearly zero!" << std::endl;
    }
}

// ============================================================================
// Enforce Boundary Conditions (Override)
// ============================================================================
void PLCProblem::enforceStrongBC(bool firstTime) {
    // First, enforce standard BCs from base class (elasticity)
    if (M_ElastPB) {
        M_ElastPB->enforceStrongBC(firstTime);
    }
    
// Then enforce displacement BC on outer wall
    if (M_useDisplacementBC) {
        enforceOuterWallDisplacementBC(firstTime);
    }
    
    // Then enforce pressure BC on outer wall (Dirichlet)
    if (M_usePolynomialBC) {
        enforceOuterWallBC(firstTime);  
    }
    
    // Add coupling matrix term
    if (firstTime && std::abs(M_gamma) > 1e-15) {
        addCouplingMatrix();
    }
    //s Add coupling matrix term if not done already
    //if (firstTime && std::abs(M_gamma) > 1e-15) {
    //    addCouplingMatrix();
    //}
}

// ============================================================================
// Enforce Polynomial Dirichlet BC on Outer Wall
// ============================================================================

void PLCProblem::enforceOuterWallBC(bool firstTime) {
    std::cout << "[PLCProblem] Enforcing polynomial BC on outer wall (region " 
              << M_outerWallRegion << ")..." << std::endl;
    
    if (!M_outerWallBCCallbackP && M_bcCoefficients.empty()) {
        std::cerr << "[PLCProblem] Warning: No BC callback or coefficients set!" << std::endl;
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
    size_type nbElastDOF = getNbElastDOF();
    size_type nbVelocityDOF = getNbVelocityDOF();
    size_type pressureOffset = nbElastDOF + nbVelocityDOF;
    
    // Get mesh and boundary region
    const getfem::mesh& mesh = mf_pressure.linked_mesh();
    
    // Check if region exists
    if (!mesh.has_region(M_outerWallRegion)) {
        std::cerr << "[PLCProblem] Warning: Region " << M_outerWallRegion 
                  << " does not exist!" << std::endl;
        return;
    }
    
    getfem::mesh_region region = mesh.region(M_outerWallRegion);
    
    // Find all pressure DOFs on the outer wall boundary
    dal::bit_vector boundary_dofs;
    
    // Iterate over faces in the outer wall region
    for (getfem::mr_visitor it(region); !it.finished(); ++it) {
        if (!it.is_face()) continue;
        
        size_type cv = it.cv();
        
        // Check if element has FEM
        getfem::pfem pf = mf_pressure.fem_of_element(cv);
        if (pf == nullptr) continue;
        
        // Get DOF indices for this element
        auto dof_indices = mf_pressure.ind_basic_dof_of_element(cv);
        
        // Mark all DOFs of this element on the boundary
        for (size_type i = 0; i < dof_indices.size(); ++i) {
            boundary_dofs.add(dof_indices[i]);
        }
    }
    
    // Apply Dirichlet BC to each boundary DOF
    size_type bc_count = 0;
    size_type ncols = gmm::mat_ncols(*matrix);
    
    for (dal::bv_visitor dof(boundary_dofs); !dof.finished(); ++dof) {
        size_type global_dof = pressureOffset + dof;
        
        // Get DOF coordinates
        bgeot::base_node pt = mf_pressure.point_of_basic_dof(dof);
        scalar_type z = pt[2];  // Assuming z is the axial direction
        
        // Evaluate polynomial BC
        scalar_type bc_value;
        if (M_outerWallBCCallbackP) {
            bc_value = M_outerWallBCCallbackP(z);
        } else {
            bc_value = evaluatePolynomial(z, M_bcCoefficients, M_z_min, M_z_max);
        }
        
        // Enforce Dirichlet BC
        if (firstTime) {
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
        scalar_type p_min = M_outerWallBCCallbackP ? M_outerWallBCCallbackP(M_z_min) 
                           : evaluatePolynomial(M_z_min, M_bcCoefficients, M_z_min, M_z_max);
        scalar_type p_max = M_outerWallBCCallbackP ? M_outerWallBCCallbackP(M_z_max) 
                           : evaluatePolynomial(M_z_max, M_bcCoefficients, M_z_min, M_z_max);
        scalar_type p_mid = M_outerWallBCCallbackP ? M_outerWallBCCallbackP((M_z_min + M_z_max) / 2.0) 
                           : evaluatePolynomial((M_z_min + M_z_max) / 2.0, M_bcCoefficients, M_z_min, M_z_max);
        
        std::cout << "  BC sample values:" << std::endl;
        std::cout << "    p(z_min=" << M_z_min << ") = " << p_min << std::endl;
        std::cout << "    p(z_mid=" << (M_z_min + M_z_max)/2.0 << ") = " << p_mid << std::endl;
        std::cout << "    p(z_max=" << M_z_max << ") = " << p_max << std::endl;
    }
}


// ============================================================================
// Print Coupling Information
// ============================================================================
void PLCProblem::printCouplingInfo() const {
    std::cout << "=====================================================" << std::endl;
    std::cout << " PLCProblem Coupling Information" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Leakage coefficient (gamma): " << M_gamma << std::endl;
    std::cout << "  Outer wall region: " << M_outerWallRegion << std::endl;
    std::cout << "  Use polynomial BC: " << (M_usePolynomialBC ? "yes" : "no") << std::endl;
    std::cout << "  Z-range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    
    if (!M_bcCoefficients.empty()) {
        std::cout << "  BC polynomial order: " << (M_bcCoefficients.size() - 1) << std::endl;
        std::cout << "  BC coefficients: [";
        for (size_t i = 0; i < M_bcCoefficients.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << M_bcCoefficients[i];
        }
        std::cout << "]" << std::endl;
    }
    
    std::cout << "  Has coupling source: " << (hasCouplingSource() ? "yes" : "no") << std::endl;
    if (hasCouplingSource()) {
        std::cout << "  Coupling source norm: " << gmm::vect_norm2(*M_couplingSource) << std::endl;
    }
    
    std::cout << "  Pressure mass matrix built: " << (M_pressureMassBuilt ? "yes" : "no") << std::endl;
    std::cout << "  Coupling matrix added: " << (M_couplingMatrixAdded ? "yes" : "no") << std::endl;
    std::cout << "=====================================================" << std::endl;
}