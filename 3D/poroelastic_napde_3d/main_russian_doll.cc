// ============================================================================
// main_russian_doll.cc - Main driver for dual-porosity (PV + PLC) simulation
// ============================================================================
/**
 * @file main_russian_doll.cc
 * @brief Main driver for Russian Doll (nested porosity) poroelasticity simulation
 * 
 * This program solves the coupled dual-porosity Biot equations with one-way coupling:
 * 
 * PV (Vascular Porosity) - Outer domain, solved first:
 *   K u_v - alpha_v C p_v = R_v
 *   alpha_v C^T u_dot_v + H_v p_v + (1/M_v) M p_dot_v = Q_v
 * 
 * PLC (Lacuno-canalicular Porosity) - Inner domain:
 *   K u_l - alpha_l C p_l = R_l
 *   alpha_l C^T u_dot_l + (H_l + gamma*M) p_l + (1/M_l) M p_dot_l = Q_l + gamma*M*p_v
 *   
 *   with BCs on outer wall:
 *   - Displacement: u_l = u_v(z)  (Dirichlet, from polynomial fit)
 *   - Pressure: p_l = p_v(z)      (Neumann in mixed form, from polynomial fit)
 * 
 * Time stepping workflow:
 * 1. Assemble PV matrix (once at start)
 * 2. For each time step:
 *    a. Clear PV RHS, assemble PV RHS, apply PV BCs
 *    b. Solve PV
 *    c. Interpolate PV solution -> polynomial coefficients
 *    d. Update PLC BCs from polynomials
 *    e. Clear PLC RHS, assemble PLC RHS (includes +gamma*M*p_v), apply PLC BCs
 *    f. Solve PLC
 *    g. Update solutions, export
 * 
 * Usage:
 *   ./main_russian_doll -f data_russian.txt
 */

#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#include <fenv.h>
#endif

// Core includes
#include "include/Core.h"
#include "include/Bulk.h"
#include "include/BC.h"
#include "include/FEM.h"
#include "include/LinearSystem.h"
#include "include/TimeLoop.h"

// Problem classes
#include "include/DarcyProblemT.h"
#include "include/ElastProblem.h"
#include "include/CoupledProblem.h"
#include "include/PLCProblem.h"
#include "include/RussianDollProblem.h"

// Utilities
#include "include/UsefulFunctions.h"

#include <sys/stat.h>

// Helper function to create directory
void createDirectory(const std::string& path) {
    mkdir(path.c_str(), 0755);
}

int main(int argc, char *argv[]) {
    
    // Enable floating point exceptions
    #ifdef GETFEM_HAVE_FEENABLEEXCEPT
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif
    
    // ========================================================================
    // Parse Command Line Arguments
    // ========================================================================
    
    GetPot command_line(argc, argv);
    const std::string data_file_name = command_line.follow("data_russian.txt", 2, "-f", "--file");
    
    std::cout << "========================================================" << std::endl;
    std::cout << "  Russian Doll Dual-Porosity Poroelasticity Simulation  " << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "Data file: " << data_file_name << std::endl;
    
    GetPot dataFile(data_file_name.data());
    
    // Output directories
    const std::string vtkFolder = "output_vtk/";
    //const std::string vtkFolder = dataFile("output/folder", "./output_vtk");
    //createDirectory(vtkFolder);
    //createDirectory(vtkFolder + "/pv");
    //createDirectory(vtkFolder + "/plc");
    
    // ========================================================================
    // Phase 1: Create Domains and Time Loop
    // ========================================================================
    
    std::cout << "\n=== Phase 1: Domain Setup ===" << std::endl;
    
    // Create bulk domains
    std::cout << "Creating PV (outer) bulk domain..." << std::endl;
    Bulk domainPV(dataFile, "bulkDataPV/");
    
    std::cout << "Creating PLC (inner) bulk domain..." << std::endl;
    Bulk domainPLC(dataFile, "bulkDataPLC/");
    
    // Export meshes
    domainPV.exportMesh(vtkFolder + "/meshPV.vtk");
    domainPLC.exportMesh(vtkFolder + "/meshPLC.vtk");
    std::cout << "Meshes exported to " << vtkFolder << std::endl;
    
    // Create time loop
    std::cout << "Setting up time discretization..." << std::endl;
    TimeLoop timeLoop(dataFile);
    std::cout << "  dt = " << timeLoop.dt() << std::endl;
    std::cout << "  Tend = " << timeLoop.Tend() << std::endl;
    std::cout << "  N steps = " << timeLoop.Nstep() << std::endl;
    
    // ========================================================================
    // Phase 2: Create Sub-Problems
    // ========================================================================
    
    std::cout << "\n=== Phase 2: Creating Sub-Problems ===" << std::endl;
    
    // --- PV Sub-Problems ---
    std::cout << "Creating PV Darcy problem..." << std::endl;
    DarcyProblemT darcyPV(dataFile, &domainPV, "bulkDataPV/");
    
    std::cout << "Creating PV elasticity problem..." << std::endl;
    ElastProblem elastPV(dataFile, &domainPV, "bulkDataPV/");
    
    // --- PLC Sub-Problems ---
    std::cout << "Creating PLC Darcy problem..." << std::endl;
    DarcyProblemT darcyPLC(dataFile, &domainPLC, "bulkDataPLC/");
    
    std::cout << "Creating PLC elasticity problem..." << std::endl;
    ElastProblem elastPLC(dataFile, &domainPLC, "bulkDataPLC/");
    
    // ========================================================================
    // Phase 3: Create Coupled Problems
    // ========================================================================
    
    std::cout << "\n=== Phase 3: Creating Coupled Problems ===" << std::endl;
    
    // PV Coupled Problem (standard CoupledProblem)
    std::cout << "Creating PV coupled problem..." << std::endl;
    CoupledProblem coupledPV(dataFile, &domainPV, &timeLoop, "bulkDataPV/");
    
    // PLC Coupled Problem (PLCProblem with polynomial BC support)
    std::cout << "Creating PLC coupled problem..." << std::endl;
    PLCProblem coupledPLC(dataFile, &domainPLC, &timeLoop, "bulkDataPLC/");
    
    // ========================================================================
    // Phase 4: Initialize Sub-Problems
    // ========================================================================
    
    std::cout << "\n=== Phase 4: Initializing Sub-Problems ===" << std::endl;
    
    // Initialize PV
    std::cout << "Initializing PV problems..." << std::endl;
    darcyPV.initialize();
    elastPV.initialize();
    
    // Initialize PLC
    std::cout << "Initializing PLC problems..." << std::endl;
    darcyPLC.initialize();
    elastPLC.initialize();
    
    // Register sub-problems with coupled systems
    std::cout << "Registering sub-problems..." << std::endl;
    coupledPV.addDarcyPB(&darcyPV);
    coupledPV.addElastPB(&elastPV);
    
    coupledPLC.addDarcyPB(&darcyPLC);
    coupledPLC.addElastPB(&elastPLC);
    
    // ========================================================================
    // Phase 5: Create Linear Systems
    // ========================================================================
    
    std::cout << "\n=== Phase 5: Creating Linear Systems ===" << std::endl;
    
    LinearSystem sysPV(dataFile, "extra/solver/");
    LinearSystem sysPLC(dataFile, "extra/solver/");
    
    // Register DOFs with systems
    std::cout << "Registering DOFs for PV system..." << std::endl;
    coupledPV.addToSys(&sysPV);
    
    std::cout << "Registering DOFs for PLC system..." << std::endl;
    coupledPLC.addToSys(&sysPLC);
    
    // ========================================================================
    // Phase 6: Create Russian Doll Manager
    // ========================================================================
    
    std::cout << "\n=== Phase 6: Creating Russian Doll Manager ===" << std::endl;
    
    RussianDollProblem russianDoll(dataFile, &domainPV, &domainPLC, &timeLoop);
    russianDoll.setPVProblem(&coupledPV);
    russianDoll.setPLCProblem(&coupledPLC);
    russianDoll.initialize();
    
    // ========================================================================
    // Phase 7: Assemble Matrices (once, they don't change)
    // ========================================================================
    
    std::cout << "\n=== Phase 7: Assembling Matrices ===" << std::endl;
    
    // Assemble PV matrix
    std::cout << "Assembling PV system matrix..." << std::endl;
    coupledPV.assembleMatrix();
    coupledPV.enforceStrongBC(true);  // Modify matrix for Dirichlet BCs
    coupledPV.addSubSystems();
    sysPV.saveMatrix("matrixPV.mm");
    
    // Assemble PLC matrix
    // Note: The coupling matrix (+gamma*M) is added in enforceStrongBC(true)
    std::cout << "Assembling PLC system matrix..." << std::endl;
    coupledPLC.assembleMatrix();
    coupledPLC.enforceStrongBC(true);  // Modify matrix + add coupling matrix
    coupledPLC.addSubSystems();
    sysPLC.saveMatrix("matrixPLC.mm");
    
    // ========================================================================
    // Phase 8: Export Initial Condition
    // ========================================================================
    
    std::cout << "\n=== Phase 8: Initial Condition (t = 0) ===" << std::endl;
    
    russianDoll.exportVtk(vtkFolder, 0);
    
    // Update for first time step
    coupledPV.updateSol();
    coupledPLC.updateSol();
    
    // ========================================================================
    // Phase 9: Time Integration Loop
    // ========================================================================
    
    std::cout << "\n=== Phase 9: Time Integration ===" << std::endl;
    std::cout << "Starting time loop for " << timeLoop.Nstep() << " time steps..." << std::endl;
    
    for (size_type tt = 0; tt < timeLoop.Nstep(); ++tt) {
        
        std::cout << "\n======================================================" << std::endl;
        std::cout << " Time Step " << (tt + 1) << " / " << timeLoop.Nstep() 
                  << " (t = " << timeLoop.time() + timeLoop.dt() << ")" << std::endl;
        std::cout << "======================================================" << std::endl;
        
        // Advance time
        timeLoop.advance();
        
        // ====================================================================
        // Step A: Solve PV Problem
        // ====================================================================
        std::cout << "\n--- Step A: Solving PV Problem ---" << std::endl;
        
        // Clear RHS
        sysPV.cleanRHS();
        coupledPV.clearSubSystemsRHS();
        
        // Assemble RHS for current time
        coupledPV.assembleRHS();
        
        // Apply BCs (RHS only, matrix already modified)
        coupledPV.enforceStrongBC(false);
        
        // Add sub-system contributions
        coupledPV.addSubSystemsRHS();
        
        // Solve
        coupledPV.solve();
        
        std::cout << "PV solution complete." << std::endl;
        
        // ====================================================================
        // Step B: Interpolate PV -> PLC
        // ====================================================================
        std::cout << "\n--- Step B: Interpolating PV -> PLC ---" << std::endl;
        
        // This extracts PV solution, fits polynomials, and updates PLC BCs
        russianDoll.interpolatePVtoPLC();
        
        // Export interpolation data for debugging
        russianDoll.exportInterpolationData(vtkFolder, tt + 1);
        
        // ====================================================================
        // Step C: Solve PLC Problem
        // ====================================================================
        std::cout << "\n--- Step C: Solving PLC Problem ---" << std::endl;
        
        // Clear RHS
        sysPLC.cleanRHS();
        coupledPLC.clearSubSystemsRHS();
        
        // Assemble RHS for current time
        // This includes the coupling term +gamma*M*p_v
        coupledPLC.assembleRHS();
        
        // Apply BCs (RHS only)
        // The polynomial BCs are already set via callbacks in BC class
        coupledPLC.enforceStrongBC(false);
        
        // Add sub-system contributions
        coupledPLC.addSubSystemsRHS();
        
        // Solve
        coupledPLC.solve();
        
        std::cout << "PLC solution complete." << std::endl;
        
        // ====================================================================
        // Step D: Post-processing
        // ====================================================================
        std::cout << "\n--- Step D: Post-processing ---" << std::endl;
        
        // Update solutions for next time step
        coupledPV.updateSol();
        coupledPLC.updateSol();
        
        // Compute errors (if exact solutions available)
        std::vector<scalar_type> errors = russianDoll.computeErrors(timeLoop.time());
        if (errors[0] >= 0) {
            std::cout << "  PV pressure error: " << errors[0] << std::endl;
        }
        if (errors[2] >= 0) {
            std::cout << "  PLC pressure error: " << errors[2] << std::endl;
        }
        
        // Export VTK
        std::cout << "Exporting solutions..." << std::endl;
        russianDoll.exportVtk(vtkFolder, tt + 1);
        
        std::cout << "Time step " << (tt + 1) << " completed." << std::endl;
    }
    
    // ========================================================================
    // Finalization
    // ========================================================================
    
    std::cout << "\n========================================================" << std::endl;
    std::cout << "  Simulation completed successfully!" << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "Results exported to: " << vtkFolder << std::endl;
    std::cout << "Total time steps: " << timeLoop.Nstep() << std::endl;
    std::cout << "Final time: " << timeLoop.time() << std::endl;
    
    // Print final coupling info
    russianDoll.printCouplingInfo();
    
    return 0;
}