// ============================================================================
// RussianDollProblem.cpp - Implementation of dual-porosity coupling
// ============================================================================

#include "../include/RussianDollProblem.h"

// ============================================================================
// Constructor
// ============================================================================
RussianDollProblem::RussianDollProblem(const GetPot& dataFile, 
                                       Bulk* bulkPV, 
                                       Bulk* bulkPLC,
                                       TimeLoop* time)
    : M_bulkPV(bulkPV),
      M_bulkPLC(bulkPLC),
      M_time(time),
      M_pvProblem(nullptr),
      M_plcProblem(nullptr),
      M_gamma(1.0e-8),
      M_pv_on_PLC(nullptr),
      M_section("russian_doll/"),
      M_initialized(false)
{
    // Read leakage coefficient
    M_gamma = dataFile((M_section + "leakage_coefficient").c_str(), 1.0e-8);
    
    // Create interpolation manager
    M_interpManager = std::make_unique<InterpolationManager>(dataFile, bulkPV, bulkPLC);
    
    std::cout << "=== RussianDollProblem Created ===" << std::endl;
    std::cout << "  Leakage coefficient (gamma): " << M_gamma << std::endl;
}

// ============================================================================
// Setup Methods
// ============================================================================
void RussianDollProblem::setPVProblem(CoupledProblem* pvProblem) {
    M_pvProblem = pvProblem;
    std::cout << "[RussianDoll] PV problem set" << std::endl;
}

void RussianDollProblem::setPLCProblem(PLCProblem* plcProblem) {
    M_plcProblem = plcProblem;
    // Pass leakage coefficient to PLC problem
    M_plcProblem->setLeakageCoefficient(M_gamma);
    std::cout << "[RussianDoll] PLC problem set" << std::endl;
}

void RussianDollProblem::initialize() {
    if (!M_pvProblem || !M_plcProblem) {
        std::cerr << "[RussianDoll] Error: Both PV and PLC problems must be set!" << std::endl;
        return;
    }
    
    // Initialize interpolation manager
    M_interpManager->initialize();
    
    // Allocate solution transfer vector
    size_type nbPressureDOF_PLC = M_plcProblem->getNbPressureDOF();
    M_pv_on_PLC.reset(new scalarVector_Type(nbPressureDOF_PLC, 0.0));  
    M_initialized = true;
    
    std::cout << "[RussianDoll] Initialized successfully" << std::endl;
    std::cout << "  PV pressure DOFs: " << M_pvProblem->getNbPressureDOF() << std::endl;
    std::cout << "  PLC pressure DOFs: " << nbPressureDOF_PLC << std::endl;
}

// ============================================================================
// Solution Methods
// ============================================================================
void RussianDollProblem::solveTimeStep() {
    if (!M_initialized) {
        std::cerr << "[RussianDoll] Error: Not initialized!" << std::endl;
        return;
    }
    
    std::cout << "\n[RussianDoll] === Solving time step ===" << std::endl;
    
    // Step 1: Solve PV problem (with p_l = 0, one-way coupling)
    solvePV();
    
    // Step 2: Interpolate PV pressure to PLC mesh
    interpolatePVtoPLC();
    
    // Step 3: Solve PLC problem with coupling source
    solvePLC();
    
    std::cout << "[RussianDoll] === Time step complete ===" << std::endl;
}

void RussianDollProblem::solvePV() {
    std::cout << "[RussianDoll] Solving PV problem..." << std::endl;
    
    // Clear and assemble PV system
    M_pvProblem->clearSubSystems();
    M_pvProblem->clearSubSystemsRHS();
    
    M_pvProblem->assembleMatrix();
    M_pvProblem->assembleRHS();
    M_pvProblem->addSubSystems();
    M_pvProblem->addSubSystemsRHS();
    M_pvProblem->enforceStrongBC(true);
    
    // Solve
    M_pvProblem->solve();
}

void RussianDollProblem::solvePLC() {
    std::cout << "[RussianDoll] Solving PLC problem..." << std::endl;
    
    // Set coupling source from PV
    M_plcProblem->setCouplingSource(M_pv_on_PLC);
    
    // Clear and assemble PLC system
    M_plcProblem->clearSubSystems();
    M_plcProblem->clearSubSystemsRHS();
    
    M_plcProblem->assembleMatrix();
    M_plcProblem->assembleRHS();  // This now includes the coupling term
    M_plcProblem->addSubSystems();
    M_plcProblem->addSubSystemsRHS();
    M_plcProblem->enforceStrongBC(true);
    
    // Solve
    M_plcProblem->solve();
}

void RussianDollProblem::interpolatePVtoPLC() {
    std::cout << "[RussianDoll] Interpolating PV -> PLC..." << std::endl;
    
    // Get PV pressure solution
    scalarVectorPtr_Type pvPressure = M_pvProblem->getPressure();
    
    if (!pvPressure || pvPressure->size() == 0) {
        std::cerr << "[RussianDoll] Warning: PV pressure not available!" << std::endl;
        return;
    }
    
    // Get mesh_fem objects for interpolation
    const getfem::mesh_fem& mf_pv = M_pvProblem->getMfPressure();
    const getfem::mesh_fem& mf_plc = M_plcProblem->getMfPressure();
    
    // Perform interpolation using the configured approach
    M_interpManager->interpolate(pvPressure, mf_pv, M_pv_on_PLC, mf_plc);
    
    std::cout << "[RussianDoll] Interpolation complete (PV norm: " 
              << gmm::vect_norm2(*pvPressure) << ", on PLC norm: "
              << gmm::vect_norm2(*M_pv_on_PLC) << ")" << std::endl;
}

void RussianDollProblem::updateSolutions() {
    std::cout << "[RussianDoll] Updating solutions..." << std::endl;
    
    M_pvProblem->updateSol();
    M_plcProblem->updateSol();
}

// ============================================================================
// Output Methods
// ============================================================================
void RussianDollProblem::exportVtk(const std::string& folder, int frame) {
    std::cout << "[RussianDoll] Exporting VTK (frame " << frame << ")..." << std::endl;
    
    // Export PV solution
    M_pvProblem->exportVtk(folder + "/pv", "all", frame);
    
    // Export PLC solution
    M_plcProblem->exportVtk(folder + "/plc", "all", frame);
}

std::vector<scalar_type> RussianDollProblem::computeErrors(scalar_type time) {
    std::vector<scalar_type> errors(4, 0.0);
    
    // Compute PV errors
    bgeot::base_node pvErrors = M_pvProblem->computeError(time);
    errors[0] = pvErrors[0];  // PV pressure error
    errors[1] = pvErrors[1];  // PV displacement error
    
    // Compute PLC errors
    bgeot::base_node plcErrors = M_plcProblem->computeError(time);
    errors[2] = plcErrors[0];  // PLC pressure error
    errors[3] = plcErrors[1];  // PLC displacement error
    
    std::cout << "[RussianDoll] Errors at t=" << time << ":" << std::endl;
    std::cout << "  PV:  p_error=" << errors[0] << ", u_error=" << errors[1] << std::endl;
    std::cout << "  PLC: p_error=" << errors[2] << ", u_error=" << errors[3] << std::endl;
    
    return errors;
}