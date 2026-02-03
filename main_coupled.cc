
// ============================================================================
// main_coupled.cc - Main driver for coupled poroelasticity simulation
// ============================================================================
/**
 * @file main_coupled.cc
 * @brief Main driver for coupled poroelasticity simulation
 * 
 * This program solves time-dependent Biot's equations (coupled Darcy flow
 * and linear elasticity) on 3D domains with optional fracture networks.
 * 
 * Physical Model:
 * - Darcy flow: K∇p = -q (permeability K, pressure p, flux q)
 * - Elasticity: σ = λtr(ε)I + 2με - αpI (stress σ, strain ε, Biot coefficient α)
 * - Coupling: ∂/∂t(ζp + α∇·u) = ∇·q (storage ζ, displacement u)
 * 
 * Numerical Method:
 * - Mixed finite elements for Darcy (RT0 for velocity, DG for pressure)
 * - Standard FEM for elasticity (continuous Lagrange elements)
 * - Backward Euler time discretization
 * - Monolithic coupling
 * - Strong enforcement of Dirichlet BC
 * 
 * Usage:
 *   ./main -f datafile.txt
 *   ./main --file datafile.txt
 */
#ifdef USE_MUMPS
    #include <mpi.h>
#endif

// Enable floating point exception handling if available
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#include <fenv.h>
#endif



// Core include - provides all GetFem++, GMM++, and standard library includes
#include "include/Core.h"

// Domain and material data classes
#include "include/Bulk.h"
#include "include/BC.h"

// Finite element wrapper
#include "include/FEM.h"

// Linear system assembly and solution
#include "include/LinearSystem.h"

// Time stepping control
#include "include/TimeLoop.h"

// Problem-specific classes
#include "include/DarcyProblemT.h"
#include "include/ElastProblem.h"
#include "include/CoupledProblem.h"

// Utility functions
#include "include/UsefulFunctions.h"

// Sparse matrix type definitions for compatibility
typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;

int main(int argc, char *argv[]) {
    // Enable floating point exception detection (helps catch NaN/Inf errors)
    #ifdef GETFEM_HAVE_FEENABLEEXCEPT
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif
    
    #ifdef USE_MUMPS
    MPI_Init(&argc, &argv);  // Needed by MUMPS library
    #endif
    // ========================================================================
    // Parse Command Line Arguments
    // ========================================================================
    
    GetPot command_line(argc, argv);
    const std::string data_file_name = command_line.follow("data", 2, "-f", "--file");
    
    std::cout << "========================================" << std::endl;
    std::cout << "  Coupled Poroelasticity Simulation    " << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Reading data file: " << data_file_name << std::endl;
    
    // Load parameter file
    GetPot dataFile(data_file_name.data());
    

    // Output directory for VTK files
    const std::string vtkFolder = "output_vtk/";
    
    // ========================================================================
    // Setup Phase: Initialize Domain, Time Loop, and Physical Problems
    // ========================================================================
    std::cout << "\n--- Initialization Phase ---" << std::endl;
    
    // Create bulk domain (loads mesh and material properties)
    std::cout << "Creating bulk domain..." << std::endl;
    Bulk myDomain(dataFile);
    myDomain.exportMesh(vtkFolder + "mesh.vtk");
    std::cout << "  Mesh exported to " << vtkFolder << "mesh.vtk" << std::endl;
    
    // Create time loop manager
    std::cout << "Setting up time discretization..." << std::endl;
    TimeLoop myTime(dataFile);
    std::cout << "  Time steps: " << myTime.Nstep() << std::endl;
    std::cout << "  Time step size: " << myTime.dt() << std::endl;
    std::cout << "  Final time: " << myTime.Tend() << std::endl;
    

    // Create individual physics problems
    std::cout << "Creating Darcy flow problem..." << std::endl;
    DarcyProblemT myDarcy(dataFile, &myDomain);
    
    std::cout << "Creating elasticity problem..." << std::endl;
    ElastProblem myElast(dataFile, &myDomain);
    
    // Create coupled problem manager
    std::cout << "Creating coupled problem..." << std::endl;
    CoupledProblem myCP(dataFile, &myDomain, &myTime);
    
    // Initialize individual problems (set FEM spaces, initial conditions)
    std::cout << "Initializing Darcy problem..." << std::endl;
    myDarcy.initialize();
    
    std::cout << "Initializing elasticity problem..." << std::endl;
    myElast.initialize();
    
    // Register sub-problems with coupled system
    myCP.addElastPB(&myElast);
    myCP.addDarcyPB(&myDarcy);
    std::cout << "  Sub-problems registered with coupled system" << std::endl;
    
    // ========================================================================
    // Assembly Phase: Build Global Monolithic System
    // ========================================================================
    
    std::cout << "\n--- Assembly Phase ---" << std::endl;
    
    LinearSystem mySys(dataFile, "solver/");
    
    std::cout << "Registering degrees of freedom..." << std::endl;
    myCP.addToSys(&mySys);
    
    std::cout << "Assembling system matrix..." << std::endl;
    myCP.assembleMatrix();
    
    std::cout << "Applying boundary conditions (first time)..." << std::endl;
    myCP.enforceStrongBC(true);  // Modify both matrix and RHS
    
    std::cout << "Combining sub-system matrices..." << std::endl;
    myCP.addSubSystems();
    
    // Save system matrix for inspection (useful for debugging)
    std::cout << "Saving system matrix to matrix.mm..." << std::endl;
    mySys.saveMatrix("matrix.mm");
    

    // ========================================================================
    // Initial Condition: Compute and Export
    // ========================================================================
    
    std::cout << "\n--- Initial Condition (t = 0) ---" << std::endl;
    
    bgeot::base_node err = myCP.computeError(0);
    std::cout << "Initial error: [pressure, displacement] = [" 
              << err[0] << ", " << err[1] << "]" << std::endl;
    
    std::cout << "Exporting initial solution..." << std::endl;
    myCP.exportVtk(vtkFolder, "all", 0);
    myCP.updateSol();
    
    // ========================================================================
    // Time Loop: Solve Transient Problem
    // ========================================================================
    
    std::cout << "\n--- Time Integration ---" << std::endl;
    std::cout << "Starting time loop for " << myTime.Nstep() << " time steps..." << std::endl;
    
    for (size_type tt = 0; tt < myTime.Nstep(); ++tt) {
        std::cout << "\n=== Time Step " << tt + 1 << " / " << myTime.Nstep() 
                  << " (t = " << myTime.time() + myTime.dt() << ") ===" << std::endl;
        
        // Advance time
        myTime.advance();
        
        // Clear previous RHS (matrix remains constant for linear problems)
        std::cout << "  Clearing RHS vectors..." << std::endl;
        mySys.cleanRHS();
        myCP.clearSubSystemsRHS();
        
        // Reassemble RHS for current time step
        std::cout << "  Assembling RHS..." << std::endl;
        myCP.assembleRHS();
        
        // Apply boundary conditions (RHS only, matrix already modified)
        std::cout << "  Applying boundary conditions..." << std::endl;
        myCP.enforceStrongBC(false);
        
        // Add sub-system contributions to global RHS
        myCP.addSubSystemsRHS();
        
        // Solve coupled system
        std::cout << "  Solving coupled system..." << std::endl;
        myCP.solve();
        
        // Update solution for next time step (old <- current)
        myCP.updateSol();
        
        // Compute error against exact solution (if available)
        std::cout << "  Computing error..." << std::endl;
        bgeot::base_node err = myCP.computeError(myTime.time());
        std::cout << "  Error: [pressure, displacement] = [" 
                  << err[0] << ", " << err[1] << "]" << std::endl;
        
        // Export results
        std::cout << "  Exporting solution..." << std::endl;
        myCP.exportVtk(vtkFolder, "all", tt + 1);
        
        std::cout << "  Time step completed." << std::endl;
    }
    
    // ========================================================================
    // Finalization
    // ========================================================================
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "  Simulation completed successfully!    " << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Results exported to: " << vtkFolder << std::endl;
    std::cout << "Total time steps: " << myTime.Nstep() << std::endl;
    std::cout << "Final time: " << myTime.time() << std::endl;
    

    #ifdef USE_MUMPS
    MPI_Finalize();
    #endif

    return 0;
}