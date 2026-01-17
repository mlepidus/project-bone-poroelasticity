// ============================================================================
// RussianDollProblem.cpp - Implementation of dual-porosity coupled problem
// ============================================================================

#include "../include/RussianDollProblem.h"
#include <iostream>
#include <cmath>

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
      M_initialized(false),
      M_section("russian_doll/")
{
    // Read coupling parameters
    M_gamma = dataFile((M_section + "gamma").c_str(), 1.0e-8);
    M_maxCouplingIter = dataFile((M_section + "max_coupling_iterations").c_str(), 10);
    M_couplingTolerance = dataFile((M_section + "coupling_tolerance").c_str(), 1.0e-6);
    
    // Read coupling approach
    std::string approach = dataFile((M_section + "interpolation/approach").c_str(), "mesh");
    if (approach == "line" || approach == "LINE") {
        M_couplingApproach = CouplingApproach::LINE_INTERPOLATION;
    } else {
        M_couplingApproach = CouplingApproach::MESH_INTERPOLATION;
    }
    
    // Create interpolation manager
    M_interpManager = std::make_unique<InterpolationManager>(dataFile, bulkPV, bulkPLC);
    
    std::cout << "====================================================" << std::endl;
    std::cout << " Russian Doll Problem Configuration" << std::endl;
    std::cout << "====================================================" << std::endl;
    std::cout << " Leakage coefficient (gamma): " << M_gamma << std::endl;
    std::cout << " Max coupling iterations: " << M_maxCouplingIter << std::endl;
    std::cout << " Coupling tolerance: " << M_couplingTolerance << std::endl;
    std::cout << " Coupling approach: " 
              << (M_couplingApproach == CouplingApproach::LINE_INTERPOLATION 
                  ? "LINE" : "MESH") << std::endl;
    std::cout << "====================================================" << std::endl;
}

// ============================================================================
// Setup Methods
// ============================================================================
void RussianDollProblem::setPVProblem(CoupledProblem* pvProblem) {
    M_pvProblem = pvProblem;
    std::cout << "[RussianDoll] PV problem registered" << std::endl;
}

void RussianDollProblem::setPLCProblem(PLCProblem* plcProblem) {
    M_plcProblem = plcProblem;
    M_plcProblem->setLeakageCoefficient(M_gamma);
    std::cout << "[RussianDoll] PLC problem registered" << std::endl;
}

void RussianDollProblem::initialize() {
    if (M_pvProblem == nullptr || M_plcProblem == nullptr) {
        std::cerr << "[RussianDoll] Error: PV or PLC problem not set!" << std::endl;
        return;
    }
    
    // Initialize interpolation structures
    M_interpManager->initialize();
    
    // Allocate transfer vectors
    size_type nb_dof_pv_p = M_pvProblem->getDarcyProblem()->getMfPressure().nb_dof();
    size_type nb_dof_plc_p = M_plcProblem->getMfPressure().nb_dof();
    
    M_pv_on_PLC = std::make_shared<scalarVector_Type>(nb_dof_plc_p, 0.0);
    M_pl_on_PV = std::make_shared<scalarVector_Type>(nb_dof_pv_p, 0.0);
    
    // Store initial pressures (for first step, p_l = 0 on PV side)
    M_pv_old = std::make_shared<scalarVector_Type>(nb_dof_pv_p, 0.0);
    M_pl_old = std::make_shared<scalarVector_Type>(nb_dof_plc_p, 0.0);
    
    M_initialized = true;
    std::cout << "[RussianDoll] Initialized successfully" << std::endl;
}

// ============================================================================
// Solution Methods - Staggered Coupling
// ============================================================================
void RussianDollProblem::solveTimeStep() {
    if (!M_initialized) {
        std::cerr << "[RussianDoll] Error: Not initialized!" << std::endl;
        return;
    }
    
    std::cout << "\n--- Russian Doll Time Step ---" << std::endl;
    
    // Staggered iteration loop
    scalar_type residual = 1.0;
    size_type iter = 0;
    
    while (residual > M_couplingTolerance && iter < M_maxCouplingIter) {
        iter++;
        std::cout << "[RussianDoll] Coupling iteration " << iter << std::endl;
        
        // Store previous PLC pressure for convergence check
        scalarVector_Type pl_prev = *M_plcProblem->getPressure();
        
        // Step 1: Solve PV with current p_l estimate
        // For first step or if one-way coupling, p_l on PV can be zero
        assembleLeakageRHS_PV(M_pl_on_PV);
        solvePV();
        
        // Step 2: Interpolate p_v to PLC mesh
        interpolatePVtoPLC();
        
        // Step 3: Solve PLC with interpolated p_v
        assembleLeakageRHS_PLC(M_pv_on_PLC);
        solvePLC();
        
        // Step 4: (Optional) Interpolate p_l back to PV for feedback
        if (M_gamma > 1e-15) {  // Only if two-way coupling
            interpolatePLCtoPV();
        }
        
        // Compute convergence: ||p_l^{k} - p_l^{k-1}|| / ||p_l^{k}||
        scalarVector_Type diff = *M_plcProblem->getPressure();
        gmm::add(gmm::scaled(pl_prev, -1.0), diff);
        
        scalar_type norm_diff = gmm::vect_norm2(diff);
        scalar_type norm_pl = gmm::vect_norm2(*M_plcProblem->getPressure());
        
        residual = (norm_pl > 1e-15) ? norm_diff / norm_pl : norm_diff;
        
        std::cout << "[RussianDoll] Coupling residual: " << residual << std::endl;
    }
    
    if (iter >= M_maxCouplingIter) {
        std::cout << "[RussianDoll] Warning: Max coupling iterations reached!" << std::endl;
    }
    
    // Update old solutions for next time step
    updateSolutions();
}

void RussianDollProblem::solvePV() {
    std::cout << "[RussianDoll] Solving PV problem..." << std::endl;
    M_pvProblem->solve();
}

void RussianDollProblem::solvePLC() {
    std::cout << "[RussianDoll] Solving PLC problem..." << std::endl;
    M_plcProblem->solve();
}

// ============================================================================
// Interpolation Methods
// ============================================================================
void RussianDollProblem::interpolatePVtoPLC() {
    std::cout << "[RussianDoll] Interpolating PV -> PLC..." << std::endl;
    
    scalarVectorPtr_Type pv_pressure = M_pvProblem->getDarcyProblem()->getPressure();
    const getfem::mesh_fem& mf_pv = M_pvProblem->getDarcyProblem()->getMfPressure();
    const getfem::mesh_fem& mf_plc = M_plcProblem->getMfPressure();
    
    M_interpManager->interpolate(pv_pressure, mf_pv, M_pv_on_PLC, mf_plc);
}

void RussianDollProblem::interpolatePLCtoPV() {
    std::cout << "[RussianDoll] Interpolating PLC -> PV..." << std::endl;
    
    // This is trickier: PLC mesh is inside PV mesh
    // For points in PV outside PLC, we need to either:
    // 1. Set to zero
    // 2. Extrapolate
    // 3. Use only the overlap region
    
    // For now, we'll set non-overlapping regions to zero
    // This is valid if γM p_l is only active where PLC exists
    
    scalarVectorPtr_Type pl_pressure = M_plcProblem->getPressure();
    const getfem::mesh_fem& mf_plc = M_plcProblem->getMfPressure();
    const getfem::mesh_fem& mf_pv = M_pvProblem->getDarcyProblem()->getMfPressure();
    
    // Initialize to zero
    gmm::clear(*M_pl_on_PV);
    
    // For mesh interpolation, use extrapolation=0 option
    // This requires custom handling since GetFEM throws for points outside
    
    // Simple approach: for each PV DOF, check if it's inside PLC domain
    // If yes, interpolate; if no, leave as zero
    
    for (size_type i = 0; i < mf_pv.nb_dof(); ++i) {
        bgeot::base_node pt = mf_pv.point_of_basic_dof(i);
        
        // Try to find this point in PLC mesh
        try {
            // Create a temporary 0D mesh at this point
            getfem::mesh pt_mesh;
            pt_mesh.add_point(pt);
            getfem::mesh_fem mf_pt(pt_mesh, 1);
            mf_pt.set_classical_finite_element(0);
            
            std::vector<scalar_type> val_at_pt(1);
            getfem::interpolation(mf_plc, mf_pt, *pl_pressure, val_at_pt);
            (*M_pl_on_PV)[i] = val_at_pt[0];
        }
        catch (...) {
            // Point is outside PLC mesh - leave as zero
            (*M_pl_on_PV)[i] = 0.0;
        }
    }
}

void RussianDollProblem::updateSolutions() {
    // Store current solutions as "old" for next time step
    gmm::copy(*M_pvProblem->getDarcyProblem()->getPressure(), *M_pv_old);
    gmm::copy(*M_plcProblem->getPressure(), *M_pl_old);
    
    // Update individual problems
    M_pvProblem->updateSolution();
    M_plcProblem->updateSolution();
}

// ============================================================================
// Assembly Methods for Leakage Terms
// ============================================================================
void RussianDollProblem::assembleLeakageRHS_PV(const scalarVectorPtr_Type& p_l_on_PV) {
    // Assemble: +γ * M * p_l
    // Where M is the mass matrix of the PV pressure space
    
    if (gmm::vect_norm2(*p_l_on_PV) < 1e-15) {
        // Skip if p_l is essentially zero (first iteration or one-way)
        return;
    }
    
    std::cout << "[RussianDoll] Assembling leakage RHS for PV..." << std::endl;
    
    // Get mass matrix from PV problem
    // Note: This assumes CoupledProblem exposes the mass matrix
    // You may need to add a getter for this
    
    // sparseMatrixPtr_Type M_mass = M_pvProblem->getDarcyProblem()->getMassMatrix();
    // scalarVectorPtr_Type leakage_term = std::make_shared<scalarVector_Type>(p_l_on_PV->size());
    // gmm::mult(*M_mass, *p_l_on_PV, *leakage_term);
    // gmm::scale(*leakage_term, M_gamma);
    // 
    // // Add to PV RHS
    // M_pvProblem->getDarcyProblem()->addToRHS(leakage_term);
}

void RussianDollProblem::assembleLeakageRHS_PLC(const scalarVectorPtr_Type& p_v_on_PLC) {
    // Set the coupling source in PLC problem
    // The PLC problem will handle: +γ * M * p_v
    M_plcProblem->setCouplingSource(p_v_on_PLC);
}

// ============================================================================
// Output Methods
// ============================================================================
void RussianDollProblem::exportVtk(const std::string& folder, int frame) {
    std::cout << "[RussianDoll] Exporting VTK frame " << frame << std::endl;
    
    // Export PV solution
    std::string pv_folder = folder + "/pv";
    M_pvProblem->exportVtk(pv_folder, frame);
    
    // Export PLC solution
    std::string plc_folder = folder + "/plc";
    M_plcProblem->exportVtk(plc_folder, frame);
    
    // Optionally export interpolated fields for visualization
    // Export p_v on PLC mesh
    if (M_pv_on_PLC && gmm::vect_norm2(*M_pv_on_PLC) > 1e-15) {
        getfem::vtk_export exp_interp(plc_folder + "/pv_interpolated_" 
                                      + std::to_string(frame) + ".vtk");
        exp_interp.exporting(M_plcProblem->getMfPressure());
        exp_interp.write_point_data(M_plcProblem->getMfPressure(), 
                                    *M_pv_on_PLC, "pv_on_plc");
    }
}

std::vector<scalar_type> RussianDollProblem::computeErrors(scalar_type time) {
    std::vector<scalar_type> errors(4, 0.0);
    
    // PV errors (if exact solution available)
    // errors[0] = p_v error
    // errors[1] = u_v error
    
    // PLC errors
    auto [p_l_err, u_l_err] = M_plcProblem->computeError(time);
    errors[2] = p_l_err;
    errors[3] = u_l_err;
    
    return errors;
}