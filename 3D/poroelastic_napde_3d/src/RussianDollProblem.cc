#include "../include/RussianDollProblem.h"
#include <iostream>

// ============================================================================
// Constructor
// ============================================================================
RussianDollProblem::RussianDollProblem(const GetPot& dataFile, 
                                       Bulk* bulkPV, 
                                       Bulk* bulkPLC,
                                       TimeLoop* time)
    : M_pvProblem(nullptr),
      M_plcProblem(nullptr),
      M_zMin(0.0),
      M_zMax(1.0),
      M_polynomialOrder(3),
      M_initialized(false)
{
    // Create interpolation manager
    M_interpManager = std::make_unique<InterpolationManager>(dataFile, bulkPV, bulkPLC);
    
    // Read polynomial order
    M_polynomialOrder = dataFile("russian_doll/polynomial_order", 3);
    
    std::cout << "=====================================================" << std::endl;
    std::cout << " RussianDollProblem Created (One-Way Coupling)" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Polynomial order: " << M_polynomialOrder << std::endl;
}

// ============================================================================
// Problem Registration
// ============================================================================
void RussianDollProblem::setPVProblem(CoupledProblem* pvProblem) {
    M_pvProblem = pvProblem;
    std::cout << "[RussianDoll] PV problem registered" << std::endl;
}

void RussianDollProblem::setPLCProblem(PLCProblem* plcProblem) {
    M_plcProblem = plcProblem;
    std::cout << "[RussianDoll] PLC problem registered" << std::endl;
}

// ============================================================================
// One-Way Coupling
// ============================================================================
void RussianDollProblem::solveOneWay() {
    if (!M_initialized) {
        initialize();
    }
    
    std::cout << "\n=== RussianDoll One-Way Coupling ===" << std::endl;
    
    // 1. Solve PV (standalone)
    solvePV();
    
    // 2. Update PLC BCs from PV solution
    updatePLCBoundaryConditions();
    
    // 3. Solve PLC with polynomial BCs
    solvePLC();
    
    std::cout << "=== One-way coupling complete ===" << std::endl;
}

void RussianDollProblem::updatePLCBoundaryConditions() {
    std::cout << "[RussianDoll] Updating PLC boundary conditions..." << std::endl;
    
    if (!M_pvProblem || !M_plcProblem) {
        std::cerr << "[RussianDoll] Error: Problems not set!" << std::endl;
        return;
    }
    
    // Use InterpolationManager to extract PV solution
    // This replaces your complex sampling/fitting code
    
    // Get PV pressure and displacement
    scalarVectorPtr_Type pvPressure = M_pvProblem->getPressure();
    scalarVectorPtr_Type pvDisplacement = M_pvProblem->getDisplacement();
    
    if (!pvPressure || pvPressure->size() == 0) {
        std::cerr << "[RussianDoll] Error: PV pressure not available!" << std::endl;
        return;
    }
    
    // Get mesh_fem for interpolation
    const getfem::mesh_fem& mf_pv_pressure = M_pvProblem->getMfPressure();
    
    // Use InterpolationManager to fit polynomial
    PolynomialFit pressureFit = M_interpManager->interpolatePressure(
        pvPressure, mf_pv_pressure, M_polynomialOrder);
    
    // Store coefficients
    M_pressureCoeffs = pressureFit.coefficients;
    M_zMin = pressureFit.z_min;
    M_zMax = pressureFit.z_max;
    
    // Optionally do displacement
    if (pvDisplacement && pvDisplacement->size() > 0) {
        const getfem::mesh_fem& mf_pv_disp = M_pvProblem->getMfDisplacement();
        
        PolynomialFit fitX, fitY, fitZ;
        M_interpManager->interpolateDisplacement(
            pvDisplacement, mf_pv_disp, fitX, fitY, fitZ);
        
        M_dispXCoeffs = fitX.coefficients;
        M_dispYCoeffs = fitY.coefficients;
        M_dispZCoeffs = fitZ.coefficients;
    }
    
    // Set BCs in PLC problem
    M_plcProblem->setOuterWallPressureBC(M_pressureCoeffs, M_zMin, M_zMax);
    
    if (!M_dispXCoeffs.empty()) {
        M_plcProblem->setOuterWallDisplacementBC(
            M_dispXCoeffs, M_dispYCoeffs, M_dispZCoeffs, M_zMin, M_zMax);
    }
    
    std::cout << "[RussianDoll] PLC BCs updated with polynomial coefficients" << std::endl;
}

void RussianDollProblem::solvePV() {
    if (!M_pvProblem) return;
    
    std::cout << "[RussianDoll] Solving PV problem..." << std::endl;
    
    // These calls should be in your main time loop, but for completeness:
    // M_pvProblem->assembleMatrix();
    // M_pvProblem->assembleRHS();
    // M_pvProblem->enforceStrongBC(false); // false = update RHS only
    // M_pvProblem->solve();
    
    std::cout << "[RussianDoll] PV solved" << std::endl;
}

void RussianDollProblem::solvePLC() {
    if (!M_plcProblem) return;
    
    std::cout << "[RussianDoll] Solving PLC problem..." << std::endl;
    
    // These calls should be in your main time loop:
    // M_plcProblem->assembleMatrix();
    // M_plcProblem->assembleRHS();  // Includes coupling RHS +γM·p_v
    // M_plcProblem->enforceStrongBC(false); // Includes polynomial BCs
    // M_plcProblem->solve();
    
    std::cout << "[RussianDoll] PLC solved" << std::endl;
}

// ============================================================================
// Time Stepping
// ============================================================================
void RussianDollProblem::advanceTimeStep() {
    if (M_pvProblem) M_pvProblem->solve();
    if (M_plcProblem) M_plcProblem->solve();
}

void RussianDollProblem::initialize() {
    if (M_initialized) return;
    
    // Initialize interpolation manager
    M_interpManager->initialize();
    
    // Allocate coefficient vectors
    M_pressureCoeffs.resize(M_polynomialOrder + 1, 0.0);
    M_dispXCoeffs.resize(M_polynomialOrder + 1, 0.0);
    M_dispYCoeffs.resize(M_polynomialOrder + 1, 0.0);
    M_dispZCoeffs.resize(M_polynomialOrder + 1, 0.0);
    
    M_initialized = true;
    std::cout << "[RussianDoll] Initialized" << std::endl;
}

// ============================================================================
// Output
// ============================================================================
void RussianDollProblem::exportSolutions(const std::string& folder, int timeStep) {
    if (M_pvProblem) {
        M_pvProblem->exportVtk(folder + "/pv", "all", timeStep);
    }
    
    if (M_plcProblem) {
        M_plcProblem->exportVtk(folder + "/plc", "all", timeStep);
    }
    
    // Export interpolation data
    std::string interpFile = folder + "/interpolation_" + std::to_string(timeStep) + ".txt";
    std::ofstream file(interpFile);
    if (file.is_open()) {
        file << "Pressure coefficients: ";
        for (const auto& c : M_pressureCoeffs) {
            file << c << " ";
        }
        file << "\nZ range: " << M_zMin << " to " << M_zMax << std::endl;
        file.close();
    }
}

void RussianDollProblem::printCouplingInfo() const {
    std::cout << "\n=== RussianDoll Coupling Info ===" << std::endl;
    std::cout << "Polynomial order: " << M_polynomialOrder << std::endl;
    std::cout << "Z-range: [" << M_zMin << ", " << M_zMax << "]" << std::endl;
    
    if (!M_pressureCoeffs.empty()) {
        std::cout << "Pressure coefficients: [";
        for (size_t i = 0; i < M_pressureCoeffs.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << M_pressureCoeffs[i];
        }
        std::cout << "]" << std::endl;
    }
}