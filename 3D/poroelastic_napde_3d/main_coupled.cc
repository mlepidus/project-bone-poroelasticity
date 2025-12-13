// ============================================================================
// Main Program
// ============================================================================

/**
 * @file main.cpp
 * @brief Main driver for coupled poroelasticity simulation
 * 
 * This program solves time-dependent Biot's equations (coupled Darcy flow
 * and linear elasticity) on 2D bulk domains with optional fracture networks.
 * 
 * Physical Model:
 * - Darcy flow: K∇p = -q (permeability K, pressure p, flux q)
 * - Elasticity: σ = λtr(ε)I + 2με - αpI (stress σ, strain ε, Biot coefficient α)
 * - Coupling: ∂/∂t(ζp + α∇·u) = ∇·q (storage ζ, displacement u)
 * 
 * Numerical Method:
 * - Mixed finite elements for Darcy (RT/BDM for velocity, DG for pressure)
 * - Standard FEM for elasticity (continuous Lagrange elements)
 * - Backward Euler time discretization
 * - Monolithic or staggered coupling
 * - Nitsche's method for essential boundary conditions
 * 
 * Usage:
 *   ./executable -f datafile.dat
 *   ./executable --file datafile.dat
 */

#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#include <fenv.h>  // Enable floating point exception handling
#endif

// Sparse matrix type definitions for easy reference
typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;

int main(int argc, char *argv[]) {
    // Parse command line arguments
    GetPot command_line(argc, argv);
    const std::string data_file_name = command_line.follow("data", 2, "-f", "--file");
    
    // Load parameter file
    GetPot dataFile(data_file_name.data());
    
    const std::string section = "";
    const std::string vtkFolder = "output_vtk/";
    
    // ==================================================================
    // Setup Phase: Initialize domain, time loop, and physical problems
    // ==================================================================
    
    Bulk myDomain(dataFile);
    myDomain.exportMesh(vtkFolder + "mesh.vtk");
    
    TimeLoop myTime(dataFile);
    
    DarcyProblemT myDarcy(dataFile, &myDomain);
    ElastProblem myElast(dataFile, &myDomain);
    
    CoupledProblem myCP(dataFile, &myDomain, &myTime);
    
    // Initialize individual problems
    myDarcy.initialize();
    myElast.initialize();
    
    // Register sub-problems with coupled system
    myCP.addElastPB(&myElast);
    myCP.addDarcyPB(&myDarcy);
    
    // ==================================================================
    // Assembly Phase: Build global monolithic system
    // ==================================================================
    
    LinearSystem mySys;
    
    myCP.addToSys(&mySys);          // Register DOFs
    myCP.assembleMatrix(&mySys);    // Assemble stiffness/mass matrices
    myCP.enforceStrongBC(true);     // Apply Dirichlet BC (modify matrix)
    myCP.addSubSystems();           // Combine sub-system matrices
    
    // Save system matrix for inspection
    mySys.saveMatrix("matrix.mm");
    
    // Compute initial error and export initial condition
    bgeot::base_node err = myCP.computeError(0);
    myCP.exportVtk(vtkFolder, "all", 0);
    myCP.updateSol();
    
    // ==================================================================
    // Time Loop: Solve transient problem
    // ==================================================================
    
    for (size_type tt = 0; tt < myTime.Nstep(); ++tt) {
        std::cout << "Time step " << tt << std::endl;
        
        // Advance time
        myTime.advance();
        
        // Clear previous RHS (matrix remains constant if problem is linear)
        mySys.cleanRHS();
        myCP.clearSubSystemsRHS();
        
        // Reassemble RHS for current time
        myCP.assembleRHS(&mySys);
        myCP.enforceStrongBC(false);  // Apply BC to RHS only
        myCP.addSubSystemsRHS();
        
        // Solve coupled system
        myCP.solve();
        
        // Update solution for next time step
        myCP.updateSol();
        
        // Compute error and export results
        bgeot::base_node err = myCP.computeError(myTime.time());
        myCP.exportVtk(vtkFolder, "all", tt + 1);
    }
    
    return 0;
}