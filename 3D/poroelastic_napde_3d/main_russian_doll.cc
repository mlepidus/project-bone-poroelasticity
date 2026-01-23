// ============================================================================
// main_russian_doll.cc - Main driver for dual-porosity (PV + PLC) simulation
// ============================================================================
/**
 * @file main_russian_doll.cc
 * @brief Main driver for Russian Doll (nested porosity) poroelasticity simulation
 * 
 * This program solves the coupled dual-porosity Biot equations:
 * - PV (Vascular Porosity): Outer domain, solved first (one-way coupling)
 * - PLC (Lacuno-canalicular Porosity): Inner domain, receives BC from PV
 * 
 * Physical Model:
 * - PV and PLC are each governed by Biot's poroelasticity equations
 * - One-way coupling: PV pressure is interpolated to PLC outer wall as Dirichlet BC
 * - p_l = p_v(z) on the outer wall of PLC, where p_v(z) is a polynomial fit
 * 
 * Coupling Strategy:
 * 1. Solve PV problem (standalone)
 * 2. Extract PV pressure along vertical line
 * 3. Fit polynomial p_v(z) to the extracted data
 * 4. Apply p_v(z) as Dirichlet BC on PLC outer wall
 * 5. Solve PLC problem
 * 
 * Usage:
 *   ./main_russian_doll -f data_russian_doll.txt
 */

// Enable floating point exception handling if available
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

// Problem-specific classes
#include "include/DarcyProblemT.h"
#include "include/ElastProblem.h"
#include "include/CoupledProblem.h"
#include "include/PLCProblem.h"
#include "include/RussianDollProblem.h"

// Utility functions
#include "include/UsefulFunctions.h"

// Type definitions
typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;

int main(int argc, char *argv[]) {
    
    // Enable floating point exception detection
    #ifdef GETFEM_HAVE_FEENABLEEXCEPT
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif
    
    // ========================================================================
    // Parse Command Line Arguments
    // ========================================================================
    
    GetPot command_line(argc, argv);
    const std::string data_file_name = command_line.follow("data_russian_doll.txt", 2, "-f", "--file");
    
    std::cout << "========================================================" << std::endl;
    std::cout << "  Russian Doll Dual-Porosity Poroelasticity Simulation  " << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "Reading data file: " << data_file_name << std::endl;
    
    GetPot dataFile(data_file_name.data());
    
    // Output directory for VTK files
    const std::string vtkFolder = dataFile("output/folder", "output_vtk/");
    
    // ========================================================================
    // Phase 1: Setup Domain and Time Loop
    // ========================================================================
    
    std::cout << "\n=== Phase 1: Initialization ===" << std::endl;
    
    // Create bulk domains for PV (outer) and PLC (inner)
    std::cout << "Creating PV bulk domain..." << std::endl;
    Bulk myDomainPV(dataFile, "bulkDataPV/");
    
    std::cout << "Creating PLC bulk domain..." << std::endl;
    Bulk myDomainPLC(dataFile, "bulkDataPLC/");
    
    // Export meshes for visualization
    myDomainPV.exportMesh(vtkFolder + "meshPV.vtk");
    myDomainPLC.exportMesh(vtkFolder + "meshPLC.vtk");
    std::cout << "Meshes exported to " << vtkFolder << std::endl;
    
    // Create time loop manager (shared between PV and PLC)
    std::cout << "Setting up time discretization..." << std::endl;
    TimeLoop myTime(dataFile);
    std::cout << "  Time steps: " << myTime.Nstep() << std::endl;
    std::cout << "  Time step size (dt): " << myTime.dt() << std::endl;
    std::cout << "  Final time: " << myTime.Tend() << std::endl;
    
    // ========================================================================
    // Phase 2: Create Physics Problems
    // ========================================================================
    
    std::cout << "\n=== Phase 2: Creating Physics Problems ===" << std::endl;
    
    // PV Problems (Vascular Porosity)
    std::cout << "Creating PV Darcy problem..." << std::endl;
    DarcyProblemT myDarcyPV(dataFile, &myDomainPV, "bulkDataPV/");
    
    std::cout << "Creating PV elasticity problem..." << std::endl;
    ElastProblem myElastPV(dataFile, &myDomainPV, "bulkDataPV/");
    
    // PLC Problems (Lacuno-canalicular Porosity)
    std::cout << "Creating PLC Darcy problem..." << std::endl;
    DarcyProblemT myDarcyPLC(dataFile, &myDomainPLC, "bulkDataPLC/");
    
    std::cout << "Creating PLC elasticity problem..." << std::endl;
    ElastProblem myElastPLC(dataFile, &myDomainPLC, "bulkDataPLC/");
    
    // ========================================================================
    // Phase 3: Create Coupled Problem Managers
    // ========================================================================
    
    std::cout << "\n=== Phase 3: Creating Coupled Problems ===" << std::endl;
    
    // PV Coupled Problem (standard CoupledProblem)
    std::cout << "Creating PV coupled problem..." << std::endl;
    CoupledProblem myPV(dataFile, &myDomainPV, &myTime, "bulkDataPV/");
    
    // PLC Coupled Problem (PLCProblem with BC callback support)
    std::cout << "Creating PLC coupled problem..." << std::endl;
    PLCProblem myPLC(dataFile, &myDomainPLC, &myTime, "bulkDataPLC/");
    
    // Russian Doll Manager (orchestrates PV → PLC coupling)
    std::cout << "Creating Russian Doll coupling manager..." << std::endl;
    RussianDollProblem myRussianDoll(dataFile, &myDomainPV, &myDomainPLC, &myTime);
    
    // ========================================================================
    // Phase 4: Initialize Sub-Problems
    // ========================================================================
    
    std::cout << "\n=== Phase 4: Initializing Sub-Problems ===" << std::endl;
    
    // Initialize PV sub-problems
    std::cout << "Initializing PV Darcy problem..." << std::endl;
    myDarcyPV.initialize();
    
    std::cout << "Initializing PV elasticity problem..." << std::endl;
    myElastPV.initialize();
    
    // Initialize PLC sub-problems
    std::cout << "Initializing PLC Darcy problem..." << std::endl;
    myDarcyPLC.initialize();
    
    std::cout << "Initializing PLC elasticity problem..." << std::endl;
    myElastPLC.initialize();
    
    // Register sub-problems with coupled systems
    std::cout << "Registering sub-problems with coupled systems..." << std::endl;
    myPV.addDarcyPB(&myDarcyPV);
    myPV.addElastPB(&myElastPV);
    
    myPLC.addDarcyPB(&myDarcyPLC);
    myPLC.addElastPB(&myElastPLC);
    
    // Register coupled problems with Russian Doll
    myRussianDoll.setPVProblem(&myPV);
    myRussianDoll.setPLCProblem(&myPLC);
    
    // ========================================================================
    // Phase 5: Create Linear Systems and Initial Assembly
    // ========================================================================
    
    std::cout << "\n=== Phase 5: Creating Linear Systems ===" << std::endl;
    
    // Create separate linear systems for PV and PLC
    LinearSystem mySysPV(dataFile, "solver/");
    LinearSystem mySysPLC(dataFile, "solver/");
    
    // Register DOFs
    std::cout << "Registering DOFs for PV system..." << std::endl;
    myPV.addToSys(&mySysPV);
    
    std::cout << "Registering DOFs for PLC system..." << std::endl;
    myPLC.addToSys(&mySysPLC);
    
    // Assemble matrices (remain constant for linear problems)
    std::cout << "Assembling PV system matrix..." << std::endl;
    myPV.assembleMatrix();
    
    std::cout << "Assembling PLC system matrix..." << std::endl;
    myPLC.assembleMatrix();
    
    // Apply initial boundary conditions to matrices
    std::cout << "Applying initial BCs to PV system..." << std::endl;
    myPV.enforceStrongBC(true);
    
    // Note: PLC outer wall BC will be applied dynamically from RussianDoll
    // Standard BCs (e.g., elasticity) are applied here
    std::cout << "Applying initial BCs to PLC system (excluding outer wall)..." << std::endl;
    // The polynomial BC for outer wall will be applied in the time loop
    
    // Combine sub-system matrices
    std::cout << "Combining sub-system matrices for PV..." << std::endl;
    myPV.addSubSystems();
    
    std::cout << "Combining sub-system matrices for PLC..." << std::endl;
    myPLC.addSubSystems();
    
    // Save matrices for debugging
    mySysPV.saveMatrix("matrixPV.mm");
    mySysPLC.saveMatrix("matrixPLC.mm");
    
    // Initialize Russian Doll (sets up interpolation and BC callbacks)
    std::cout << "Initializing Russian Doll coupling..." << std::endl;
    myRussianDoll.initialize();
    
    // ========================================================================
    // Phase 6: Initial Condition
    // ========================================================================
    
    std::cout << "\n=== Phase 6: Initial Condition (t = 0) ===" << std::endl;
    
    // Export initial state
    myRussianDoll.exportVtk(vtkFolder, 0);
    
    // Update solutions for first time step
    myPV.updateSol();
    myPLC.updateSol();
    
    // ========================================================================
    // Phase 7: Time Integration Loop
    // ========================================================================
    
    std::cout << "\n=== Phase 7: Time Integration ===" << std::endl;
    std::cout << "Starting time loop for " << myTime.Nstep() << " time steps..." << std::endl;
    
    for (size_type tt = 0; tt < myTime.Nstep(); ++tt) {
        
        std::cout << "\n======================================================" << std::endl;
        std::cout << " Time Step " << tt + 1 << " / " << myTime.Nstep() 
                  << " (t = " << myTime.time() + myTime.dt() << ")" << std::endl;
        std::cout << "======================================================" << std::endl;
        
        // Advance time
        myTime.advance();
        
        // --------------------------------------------------------------------
        // Step 1: Solve PV Problem
        // --------------------------------------------------------------------
        std::cout << "\n--- Solving PV Problem ---" << std::endl;
        
        // Clear previous RHS
        mySysPV.cleanRHS();
        myPV.clearSubSystemsRHS();
        
        // Assemble RHS for current time
        myPV.assembleRHS();
        
        // Apply BCs (RHS only, matrix already modified)
        myPV.enforceStrongBC(false);
        
        // Add sub-system contributions
        myPV.addSubSystemsRHS();
        
        // Solve
        myPV.solve();
        
        std::cout << "PV solution complete." << std::endl;
        
        // --------------------------------------------------------------------
        // Step 2: Interpolate PV → PLC
        // --------------------------------------------------------------------
        std::cout << "\n--- Interpolating PV Pressure to PLC Boundary ---" << std::endl;
        
        // Extract PV pressure and fit polynomial
        myRussianDoll.interpolatePVtoPLC();
        
        // Export interpolation data for debugging
        myRussianDoll.exportInterpolationData(vtkFolder, tt + 1);
        
        // --------------------------------------------------------------------
        // Step 3: Solve PLC Problem with Polynomial BC
        // --------------------------------------------------------------------
        std::cout << "\n--- Solving PLC Problem ---" << std::endl;
        
        // Clear previous RHS
        mySysPLC.cleanRHS();
        myPLC.clearSubSystemsRHS();
        
        // Assemble RHS for current time
        myPLC.assembleRHS();
        
        // Apply BCs including polynomial BC on outer wall
        // The polynomial coefficients were set by myRussianDoll.interpolatePVtoPLC()
        myPLC.enforceStrongBC(false);
        
        // Add sub-system contributions
        myPLC.addSubSystemsRHS();
        
        // Solve
        myPLC.solve();
        
        std::cout << "PLC solution complete." << std::endl;
        
        // --------------------------------------------------------------------
        // Step 4: Update and Export
        // --------------------------------------------------------------------
        std::cout << "\n--- Post-processing ---" << std::endl;
        
        // Update solutions for next time step
        myPV.updateSol();
        myPLC.updateSol();
        
        // Compute errors (if exact solution available)
        std::cout << "Computing errors..." << std::endl;
        std::vector<scalar_type> errors = myRussianDoll.computeErrors(myTime.time());
        
        // Export VTK
        std::cout << "Exporting solutions..." << std::endl;
        myRussianDoll.exportVtk(vtkFolder, tt + 1);
        
        std::cout << "Time step " << tt + 1 << " completed." << std::endl;
    }
    
    // ========================================================================
    // Finalization
    // ========================================================================
    
    std::cout << "\n========================================================" << std::endl;
    std::cout << "  Simulation completed successfully!" << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "Results exported to: " << vtkFolder << std::endl;
    std::cout << "Total time steps: " << myTime.Nstep() << std::endl;
    std::cout << "Final time: " << myTime.time() << std::endl;
    
    return 0;
}