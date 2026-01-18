// ============================================================================
// PLCProblem.cpp - Implementation of lacuno-canalicular porosity problem
// ============================================================================

#include "../include/PLCProblem.h"

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
      M_pressureMassBuilt(false)
{
    // Read leakage coefficient from data file
    // Look in section (e.g., "plc/") first, then global
    M_gamma = dataFile((section + "leakage_coefficient").c_str(), 
                       dataFile("leakage_coefficient", 1.0e-8));
    
    std::cout << "=== PLCProblem Created ===" << std::endl;
    std::cout << "  Section: " << section << std::endl;
    std::cout << "  Leakage coefficient (gamma): " << M_gamma << std::endl;
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
    M_couplingSource = std::make_shared<scalarVector_Type>(*pv_on_PLC);
    
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
    M_pressureMass = std::make_shared<sparseMatrix_Type>(nbPressureDOF, nbPressureDOF);
    gmm::clear(*M_pressureMass);
    
    // Get pressure FEM and integration method
    FEM* pressureFEM = M_DarcyPB->getFEM("Pressure");
    
    // Build mass matrix: M_ij = ∫ N_i * N_j dV
    // Using the UsefulFunctions from your codebase
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
    
    // Also build the pressure mass matrix for coupling term
    buildPressureMassMatrix();
}

// ============================================================================
// Assemble RHS (Override)
// ============================================================================
void PLCProblem::assembleRHS() {
    // Call base class to assemble standard RHS
    CoupledProblem::assembleRHS();
    
    // Add the inter-porosity coupling term
    assembleCouplingRHS();
}

// ============================================================================
// Assemble Coupling RHS Term
// ============================================================================
void PLCProblem::assembleCouplingRHS() {
    // Add coupling term: +γ * M * p_v to the pressure equation RHS
    // This corresponds to: Q_l^N + γM p_v^N in the governing equations
    
    // Skip if no coupling source or zero leakage
    if (!hasCouplingSource()) {
        return;
    }
    
    if (std::abs(M_gamma) < 1.0e-15) {
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
    scalarVectorPtr_Type couplingRHS = std::make_shared<scalarVector_Type>(nbPressureDOF, 0.0);
    
    // Matrix-vector multiplication: couplingRHS = M * p_v
    gmm::mult(*M_pressureMass, *M_couplingSource, *couplingRHS);
    
    // Scale by leakage coefficient: couplingRHS = γ * M * p_v
    gmm::scale(*couplingRHS, M_gamma);
    
    // Add to global RHS at the correct position
    // In the coupled system, DOFs are ordered as: [Elasticity | Velocity | Pressure]
    // The pressure DOFs start at offset: nbElastDOF + nbVelocityDOF
    size_type pressureOffset = nbElastDOF + nbVelocityDOF;
    
    // Get global RHS and add coupling contribution
    scalarVectorPtr_Type globalRHS = M_Sys->getRHS();
    for (size_type i = 0; i < nbPressureDOF; ++i) {
        (*globalRHS)[pressureOffset + i] += (*couplingRHS)[i];
    }
    
    std::cout << "[PLCProblem] Coupling RHS assembled (norm = " 
              << gmm::vect_norm2(*couplingRHS) << ")" << std::endl;
}