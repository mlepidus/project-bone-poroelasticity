// ============================================================================
// PLCProblem.cc - Implementation of lacuno-canalicular porosity problem
// ============================================================================
// Simplified version that delegates polynomial BC handling to the BC class.
// 
// Workflow:
// 1. RussianDollProblem calls setOuterWallPressureBC() and setOuterWallDisplacementBC()
// 2. These methods set callbacks in the BC class (for Darcy and Elasticity)
// 3. Standard assembly and BC enforcement use the callbacks automatically
// 4. Coupling RHS (+gamma * M * p_v) is added in assembleRHS()
// 5. Coupling matrix (+gamma * M) is added in enforceStrongBC(firstTime=true)
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
      M_couplingPressure(nullptr),
      M_pressureMass(nullptr),
      M_pressureMassBuilt(false),
      M_couplingMatrixAdded(false),
      M_z_min(0.0),
      M_z_max(1.0),
      M_outerWallRegion(0)
{
    // Read leakage coefficient from data file
    M_gamma = dataFile((section + "darcy/leakage").c_str(), 
                       dataFile("leakage_coefficient", 0.0));
    
    // Read outer wall region ID
    M_outerWallRegion = dataFile((section + "outer_wall_region").c_str(), 0);
    
    std::cout << "=====================================================" << std::endl;
    std::cout << " PLCProblem Created (Simplified BC Version)" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Section: " << section << std::endl;
    std::cout << "  Leakage coefficient (gamma): " << M_gamma << std::endl;
    std::cout << "  Outer wall region: " << M_outerWallRegion << std::endl;
    std::cout << "=====================================================" << std::endl;
}

// ============================================================================
// Polynomial BC Setup (delegates to BC class)
// ============================================================================

void PLCProblem::setOuterWallPressureBC(const std::vector<scalar_type>& coefficients,
                                        scalar_type z_min, scalar_type z_max) {
    if (!M_DarcyPB) {
        std::cerr << "[PLCProblem] Error: Darcy problem not set!" << std::endl;
        return;
    }
    
    // Store for debugging/export
    M_pressureCoefficients = coefficients;
    M_z_min = z_min;
    M_z_max = z_max;
    
    // Get the BC object from Darcy problem and set the polynomial callback
    // The Darcy BC class will use this callback in BCNeum() for the outer wall region
    BC* darcyBC = M_DarcyPB->getBC();
    if (darcyBC) {
        darcyBC->setPressurePolynomial(M_outerWallRegion, coefficients, z_min, z_max);
    }
    
    std::cout << "[PLCProblem] Outer wall pressure BC set via BC class" << std::endl;
    std::cout << "  Region: " << M_outerWallRegion << std::endl;
    std::cout << "  Z range: [" << z_min << ", " << z_max << "]" << std::endl;
    std::cout << "  Polynomial order: " << (coefficients.size() - 1) << std::endl;
    
    // Also set up the coupling pressure from the same polynomial
    setCouplingPressureFromPolynomial(coefficients, z_min, z_max);
}

void PLCProblem::setOuterWallDisplacementBC(const std::vector<scalar_type>& coeffs_x,
                                            const std::vector<scalar_type>& coeffs_y,
                                            const std::vector<scalar_type>& coeffs_z,
                                            scalar_type z_min, scalar_type z_max) {
    if (!M_ElastPB) {
        std::cerr << "[PLCProblem] Error: Elasticity problem not set!" << std::endl;
        return;
    }
    
    // Store for debugging/export
    M_dispXCoefficients = coeffs_x;
    M_dispYCoefficients = coeffs_y;
    M_dispZCoefficients = coeffs_z;
    M_z_min = z_min;
    M_z_max = z_max;
    
    // Get the BC object from Elasticity problem and set the polynomial callback
    // The Elasticity BC class will use this callback in BCDiriVec() for the outer wall region
    BC* elastBC = M_ElastPB->getBC();
    if (elastBC) {
        elastBC->setDisplacementPolynomial(M_outerWallRegion, 
                                           coeffs_x, coeffs_y, coeffs_z,
                                           z_min, z_max);
    }
    
    std::cout << "[PLCProblem] Outer wall displacement BC set via BC class" << std::endl;
    std::cout << "  Region: " << M_outerWallRegion << std::endl;
    std::cout << "  Z range: [" << z_min << ", " << z_max << "]" << std::endl;
}

void PLCProblem::clearOuterWallPolynomialBCs() {
    // Clear pressure callback
    if (M_DarcyPB) {
        BC* darcyBC = M_DarcyPB->getBC();
        if (darcyBC) {
            darcyBC->clearPressureCallback(M_outerWallRegion);
        }
    }
    
    // Clear displacement callback
    if (M_ElastPB) {
        BC* elastBC = M_ElastPB->getBC();
        if (elastBC) {
            elastBC->clearDisplacementCallback(M_outerWallRegion);
        }
    }
    
    // Clear stored coefficients
    M_pressureCoefficients.clear();
    M_dispXCoefficients.clear();
    M_dispYCoefficients.clear();
    M_dispZCoefficients.clear();
    
    std::cout << "[PLCProblem] Polynomial BCs cleared for region " << M_outerWallRegion << std::endl;
}

// ============================================================================
// Coupling Methods
// ============================================================================

void PLCProblem::setLeakageCoefficient(scalar_type gamma) {
    M_gamma = gamma;
    std::cout << "[PLCProblem] Leakage coefficient set to: " << M_gamma << std::endl;
}

void PLCProblem::setCouplingPressure(const scalarVectorPtr_Type& pv_on_PLC) {
    if (!pv_on_PLC) {
        std::cerr << "[PLCProblem] Warning: Null coupling pressure provided!" << std::endl;
        return;
    }
    
    size_type expectedSize = getNbPressureDOF();
    if (pv_on_PLC->size() != expectedSize) {
        std::cerr << "[PLCProblem] Warning: Coupling pressure size mismatch! "
                  << "Expected " << expectedSize 
                  << ", got " << pv_on_PLC->size() << std::endl;
    }
    
    // Copy the coupling pressure
    M_couplingPressure.reset(new scalarVector_Type(*pv_on_PLC));
    
    std::cout << "[PLCProblem] Coupling pressure set (norm = " 
              << gmm::vect_norm2(*M_couplingPressure) << ")" << std::endl;
}

void PLCProblem::setCouplingPressureFromPolynomial(const std::vector<scalar_type>& coefficients,
                                                   scalar_type z_min, scalar_type z_max) {
    if (!M_DarcyPB) {
        std::cerr << "[PLCProblem] Error: Darcy problem not set!" << std::endl;
        return;
    }
    
    size_type nbPressureDOF = getNbPressureDOF();
    
    // Allocate coupling pressure vector
    M_couplingPressure.reset(new scalarVector_Type(nbPressureDOF, 0.0));
    
    // Get pressure FEM
    const getfem::mesh_fem& mf_pressure = *(M_DarcyPB->getFEM("Pressure")->getFEM());
    
    // Helper to evaluate polynomial
    auto evaluatePoly = [&coefficients, z_min, z_max](scalar_type z) -> scalar_type {
        if (coefficients.empty()) return 0.0;
        
        scalar_type t = 0.0;
        if (std::abs(z_max - z_min) > 1e-15) {
            t = (z - z_min) / (z_max - z_min);
        }
        t = std::max(0.0, std::min(1.0, t));
        
        scalar_type result = 0.0;
        scalar_type t_power = 1.0;
        for (size_type i = 0; i < coefficients.size(); ++i) {
            result += coefficients[i] * t_power;
            t_power *= t;
        }
        return result;
    };
    
    // Evaluate polynomial at each pressure DOF location
    for (size_type i = 0; i < nbPressureDOF; ++i) {
        bgeot::base_node pt = mf_pressure.point_of_basic_dof(i);
        scalar_type z = pt[2];  // Assuming z is axial direction
        (*M_couplingPressure)[i] = evaluatePoly(z);
    }
    
    std::cout << "[PLCProblem] Coupling pressure computed from polynomial (norm = " 
              << gmm::vect_norm2(*M_couplingPressure) << ")" << std::endl;
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
    
    // Build mass matrix using massL2Standard from DarcyOperatorsBulk
    massL2Standard(M_pressureMass, 
                   *pressureFEM,
                   *pressureFEM,
                   M_intMethod);
    
    M_pressureMassBuilt = true;
    
    std::cout << "[PLCProblem] Pressure mass matrix built (" 
              << nbPressureDOF << " x " << nbPressureDOF << ")" << std::endl;
}

// ============================================================================
// Add Coupling Matrix (+gamma * M to pressure block)
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
    // This includes:
    // - Darcy source terms
    // - Darcy Neumann BCs (including outer wall pressure via callback!)
    // - Elasticity body forces
    // - Elasticity Neumann BCs (traction)
    CoupledProblem::assembleRHS();
    
    // Add the inter-porosity coupling term: +gamma * M * p_v
    if (std::abs(M_gamma) > 1e-15 && hasCouplingPressure()) {
        assembleCouplingRHS();
    }
}

// ============================================================================
// Assemble Coupling RHS Term
// ============================================================================

void PLCProblem::assembleCouplingRHS() {
    if (!hasCouplingPressure()) {
        std::cout << "[PLCProblem] No coupling pressure - skipping coupling RHS" << std::endl;
        return;
    }
    
    if (std::abs(M_gamma) < 1e-15) {
        return;
    }
    
    if (!M_pressureMassBuilt) {
        std::cerr << "[PLCProblem] Warning: Pressure mass matrix not built! Building now..." << std::endl;
        buildPressureMassMatrix();
    }
    
    size_type nbPressureDOF = getNbPressureDOF();
    size_type nbVelocityDOF = getNbVelocityDOF();
    size_type nbElastDOF = getNbElastDOF();
    
    // Compute coupling RHS: gamma * M * p_v
    scalarVectorPtr_Type couplingRHS;
    couplingRHS.reset(new scalarVector_Type(nbPressureDOF, 0.0));
    
    // Matrix-vector multiplication: couplingRHS = M * p_v
    gmm::mult(*M_pressureMass, *M_couplingPressure, *couplingRHS);
    
    // Scale by leakage coefficient: couplingRHS = gamma * M * p_v
    gmm::scale(*couplingRHS, M_gamma);
    
    // Add to global RHS at the correct position
    // DOFs are ordered as: [Elasticity | Velocity | Pressure]
    size_type pressureOffset = nbElastDOF + nbVelocityDOF;
    
    scalarVectorPtr_Type globalRHS = M_Sys->getRHS();
    for (size_type i = 0; i < nbPressureDOF; ++i) {
        (*globalRHS)[pressureOffset + i] += (*couplingRHS)[i];
    }
    
    // Diagnostics
    scalar_type rhs_norm = gmm::vect_norm2(*couplingRHS);
    std::cout << "[PLCProblem] Coupling RHS assembled:" << std::endl;
    std::cout << "  gamma = " << M_gamma << std::endl;
    std::cout << "  ||gamma * M * p_v|| = " << rhs_norm << std::endl;
}

// ============================================================================
// Enforce Boundary Conditions (Override)
// ============================================================================

void PLCProblem::enforceStrongBC(bool firstTime) {
    // Call base class which handles:
    // - Elasticity Dirichlet BCs via Nitsche or strong enforcement
    //   (including outer wall displacement via BCDiriVec callback!)
    // - The BC class callbacks are automatically used
    CoupledProblem::enforceStrongBC(firstTime);
    
    // Add coupling matrix term (+gamma * M) on first call
    if (firstTime && std::abs(M_gamma) > 1e-15) {
        addCouplingMatrix();
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
    std::cout << "  Z-range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    
    // Check BC callbacks
    if (M_DarcyPB && M_DarcyPB->getBC()) {
        bool hasPressureCB = M_DarcyPB->getBC()->hasPressureCallback(M_outerWallRegion);
        std::cout << "  Darcy pressure callback on region " << M_outerWallRegion 
                  << ": " << (hasPressureCB ? "YES" : "NO") << std::endl;
    }
    
    if (M_ElastPB && M_ElastPB->getBC()) {
        bool hasDispCB = M_ElastPB->getBC()->hasDisplacementCallback(M_outerWallRegion);
        std::cout << "  Elasticity displacement callback on region " << M_outerWallRegion 
                  << ": " << (hasDispCB ? "YES" : "NO") << std::endl;
    }
    
    if (!M_pressureCoefficients.empty()) {
        std::cout << "  Pressure polynomial order: " << (M_pressureCoefficients.size() - 1) << std::endl;
        std::cout << "  Pressure coefficients: [";
        for (size_t i = 0; i < M_pressureCoefficients.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << M_pressureCoefficients[i];
        }
        std::cout << "]" << std::endl;
    }
    
    std::cout << "  Has coupling pressure: " << (hasCouplingPressure() ? "YES" : "NO") << std::endl;
    if (hasCouplingPressure()) {
        std::cout << "  Coupling pressure norm: " << gmm::vect_norm2(*M_couplingPressure) << std::endl;
    }
    
    std::cout << "  Pressure mass matrix built: " << (M_pressureMassBuilt ? "YES" : "NO") << std::endl;
    std::cout << "  Coupling matrix added: " << (M_couplingMatrixAdded ? "YES" : "NO") << std::endl;
    std::cout << "=====================================================" << std::endl;
}