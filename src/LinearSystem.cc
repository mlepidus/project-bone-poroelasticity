#include "../include/LinearSystem.h"

#ifdef USE_MUMPS
    #include <gmm/gmm_MUMPS_interface.h>
#endif

// ============================================================================
// Constructors
// ============================================================================

LinearSystem::LinearSystem() :
    M_Matrix(),
    M_RHS(),
    M_gotInverse(false),
    M_ndof(0),
    M_iterations(0),
    M_residual(0.0),
    M_converged(false),
    M_solverType(SolverType::SUPERLU),
    M_precondType(PreconditionerType::NONE),
    M_maxIter(1000),
    M_tolerance(1e-8),
    M_verbose(false),
    M_configured(false)  
{
    #ifdef VERBOSE
        std::cout << "LinearSystem created (unconfigured - call configureSolver() or use other constructor)" << std::endl;
        M_verbose = true;
    #endif
}

LinearSystem::LinearSystem(const GetPot& dataFile, const std::string& section) :
    M_Matrix(),
    M_RHS(),
    M_gotInverse(false),
    M_ndof(0),
    M_iterations(0),
    M_residual(0.0),
    M_converged(false),
    M_solverType(SolverType::SUPERLU),
    M_precondType(PreconditionerType::NONE),
    M_maxIter(1000),
    M_tolerance(1e-8),
    M_verbose(false),
    M_configured(false)
{
    // Automatically configure solver from data file
    configureSolver(dataFile, section);
    #ifdef VERBOSE
        std::cout << "LinearSystem created (unconfigured - call configureSolver() or use other constructor)" << std::endl;
        M_verbose = true;
    #endif

}

void LinearSystem::addToMatrix(int ndof)
{
	M_ndof=ndof;
	if (M_RHS!=NULL)
	{
		int provv_size(M_RHS->size());
		M_RHS->resize(provv_size+ndof);
		M_Sol->resize(provv_size+ndof);
		
		M_Matrix->resize(provv_size+ndof,provv_size+ndof);
	}
	else
	{
		M_RHS.reset(new scalarVector_Type (ndof));
		M_Sol.reset(new scalarVector_Type (ndof));		
		M_Matrix.reset(new sparseMatrix_Type (ndof,ndof));
	}
}



void LinearSystem::copySubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column, scalar_type scale, bool transpose)
{
	if (transpose)
	{
	gmm::copy ( gmm::transposed(gmm::scaled(*M,scale)), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (first_row,gmm::mat_ncols(*M)),
                        gmm::sub_interval (first_column, gmm::mat_nrows(*M)) ) );
	}
	else
	{
	gmm::copy ( gmm::scaled(*M,scale), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (first_row,gmm::mat_nrows(*M)),
                        gmm::sub_interval (first_column, gmm::mat_ncols(*M)) ) );
	}
}

void LinearSystem::addSubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column, scalar_type scale, bool transpose)
{
	if (transpose)
	{
	gmm::add ( gmm::transposed(gmm::scaled(*M,scale)), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (first_row,gmm::mat_ncols(*M)),
                        gmm::sub_interval (first_column, gmm::mat_nrows(*M)) ) );
	}
	else
	{
	gmm::add ( gmm::scaled(*M,scale), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (first_row,gmm::mat_nrows(*M)),
                        gmm::sub_interval (first_column, gmm::mat_ncols(*M)) ) );
	}
}

void LinearSystem::copySubVector(scalarVectorPtr_Type M, int first_row, scalar_type scale)
{
	
	gmm::copy ( gmm::scaled(*M,scale), gmm::sub_vector (*M_RHS,
                        gmm::sub_interval (first_row,M->size()) ));
}

void LinearSystem::addSubVector(scalarVectorPtr_Type M, int first_row, scalar_type scale)
{
	gmm::add ( gmm::scaled(*M,scale), gmm::sub_vector (*M_RHS,
                        gmm::sub_interval (first_row, M->size()) ));
}

void LinearSystem::extractSubVector(scalarVectorPtr_Type M, int first_row, std::string where)
{
	if (where=="sol")
	{
		gmm::copy ( gmm::sub_vector (*M_Sol, gmm::sub_interval (first_row, M->size()) ), *M);
	}
	else
	{
		gmm::copy ( gmm::sub_vector (*M_RHS, gmm::sub_interval (first_row,M->size()) ), *M);
	}
}

void LinearSystem::addSubSystem(LinearSystem* small, size_type shiftRows, size_type shiftColumns)
{
	gmm::add ( *(small->getMatrix()), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (shiftRows,gmm::mat_nrows(*(small->getMatrix()))),
                        gmm::sub_interval (shiftColumns, gmm::mat_ncols(*(small->getMatrix()))) ) );
                        
        gmm::add ( *(small->getRHS()), gmm::sub_vector (*M_RHS,
                        gmm::sub_interval (shiftRows,gmm::mat_nrows(*(small->getMatrix()))) ) );
                        
}

void LinearSystem::addSubSystemRHS(LinearSystem* small, size_type shiftRows)
{
                        
        gmm::add ( *(small->getRHS()), gmm::sub_vector (*M_RHS,
                        gmm::sub_interval (shiftRows,gmm::mat_nrows(*(small->getMatrix()))) ) );
                        
}

void LinearSystem::multAddToRHS(scalarVectorPtr_Type V, int first_row, int first_column, size_type nrows, size_type ncols)
{
	int length((*V).size());	
	if (ncols==(*V).size()){
	gmm::mult_add(gmm::sub_matrix(*M_Matrix, gmm::sub_interval (first_row,nrows),gmm::sub_interval (first_column,length)),*V, gmm::sub_vector(*M_RHS,gmm::sub_interval(first_row,nrows))); 
	}
	else
	{
		std::cout << "dimension mismatch"<<std::endl;
	}
}

void LinearSystem::multAddToRHS(sparseMatrixPtr_Type M, scalarVectorPtr_Type V,  int first_rowVector, int first_rowRHS, scalar_type scale,bool transposed)
{
	if (transposed)
	{
        gmm::mult_add(gmm::scaled(gmm::transposed(*M),scale), gmm::sub_vector(*V,gmm::sub_interval(first_rowVector,gmm::mat_nrows(*M))), gmm::sub_vector(*M_RHS,gmm::sub_interval		(first_rowRHS,gmm::mat_ncols(*M)))); 
	}
	else
	{
	gmm::mult_add(gmm::scaled(*M,scale), gmm::sub_vector(*V,gmm::sub_interval(first_rowVector,gmm::mat_ncols(*M))), gmm::sub_vector(*M_RHS,gmm::sub_interval(first_rowRHS,gmm::mat_nrows(*M)))); 
	}
}

void LinearSystem::multAddToRHS(sparseMatrixPtr_Type M, scalarVector_Type& V,  int first_rowVector, int first_rowRHS, scalar_type scale,bool transposed)
{
	if (transposed)
	{
	gmm::mult_add(gmm::scaled(gmm::transposed(*M),scale),gmm::sub_vector(V,gmm::sub_interval(first_rowVector,gmm::mat_nrows(*M))), gmm::sub_vector(*M_RHS,gmm::sub_interval(first_rowRHS,gmm::mat_ncols(*M)))); 
	}
	else
	{
	gmm::mult_add(gmm::scaled(*M,scale),gmm::sub_vector(V,gmm::sub_interval(first_rowVector,gmm::mat_ncols(*M))), gmm::sub_vector(*M_RHS,gmm::sub_interval(first_rowRHS,gmm::mat_nrows(*M)))); 
	}
}


void LinearSystem::configureSolver(const GetPot& dataFile, const std::string& section)
{
    if (M_verbose)
        std::cout << "\n=== Configuring Linear Solver ===" << std::endl;
    
    // Read solver type
    std::string solverStr = dataFile((section + "type").data(), "SUPERLU");
    if (M_verbose)
    std::cout << "Solver type: " << solverStr << std::endl;
    
    if (solverStr == "SUPERLU" || solverStr == "superlu")
        M_solverType = SolverType::SUPERLU;
    else if (solverStr == "GMRES" || solverStr == "gmres")
        M_solverType = SolverType::GMRES;
    else if (solverStr == "BICGSTAB" || solverStr == "bicgstab" || solverStr == "BiCGSTAB")
        M_solverType = SolverType::BICGSTAB;
    #ifdef USE_MUMPS
    else if (solverStr == "MUMPS" || solverStr == "mumps")
        M_solverType = SolverType::MUMPS;
    #endif
    else {
        std::cerr << "WARNING: Unknown solver type '" << solverStr 
                  << "', defaulting to SUPERLU" << std::endl;
        M_solverType = SolverType::SUPERLU;
    }
    
    // Read preconditioner type
    std::string precondStr = dataFile((section + "preconditioner").data(), "NONE");
    if (M_verbose)
    std::cout << "Preconditioner: " << precondStr << std::endl;
    
    if (precondStr == "NONE" || precondStr == "none")
        M_precondType = PreconditionerType::NONE;
    else if (precondStr == "DIAGONAL" || precondStr == "diagonal" || precondStr == "jacobi")
        M_precondType = PreconditionerType::DIAGONAL;
    else if (precondStr == "ILU" || precondStr == "ilu")
        M_precondType = PreconditionerType::ILU;
    else if (precondStr == "ILUT" || precondStr == "ilut")
        M_precondType = PreconditionerType::ILUT;
    else {
        std::cerr << "WARNING: Unknown preconditioner '" << precondStr 
                  << "', defaulting to NONE" << std::endl;
        M_precondType = PreconditionerType::NONE;
    }
    
    // Read iteration parameters
    M_maxIter = dataFile((section + "maxIterations").data(), 1000);
    M_tolerance = dataFile((section + "tolerance").data(), 1e-8);
    M_verbose = dataFile((section + "verbose").data(), 1) != 0;
    if (M_verbose){
    std::cout << "Max iterations: " << M_maxIter << std::endl;
    std::cout << "Tolerance: " << M_tolerance << std::endl;
    std::cout << "Verbose: " << (M_verbose ? "true" : "false") << std::endl;
      
    }

    M_configured = true;  // Mark as configured
    std::cout << "=== Solver Configuration Complete ===\n" << std::endl;
}

// ============================================================================
// Main Solve - Uses Stored Configuration
// ============================================================================

void LinearSystem::solve()
{
    // Check if solver was configured
    if (!M_configured)
    {
        std::cerr << "WARNING: Solver not configured! Using default settings (SuperLU)." << std::endl;
        std::cerr << "Call configureSolver() or use the constructor with dataFile parameter." << std::endl;
    }
    
    M_iterations = 0;
    M_residual = 0.0;
    M_converged = false;
    
    if (M_verbose)
    {
        std::cout << "\n=== Linear System Solve ===" << std::endl;
        std::cout << "System size: " << M_ndof << " DOFs" << std::endl;
        std::cout << "Matrix non-zeros: " << gmm::nnz(*M_Matrix) << std::endl;
    }
    
    switch (M_solverType)
    {
        case SolverType::SUPERLU:
            if (M_verbose) std::cout << "Solver: SuperLU (direct)" << std::endl;
            solveDirect_SuperLU();
            break;
        
        #ifdef USE_MUMPS
        case SolverType::MUMPS:
            if (M_verbose) std::cout << "Solver: MUMPS (direct)" << std::endl;
            solveDirect_MUMPS();
            break;
        #endif
        case SolverType::GMRES:
            if (M_verbose)
            {
                std::cout << "Solver: GMRES (iterative)" << std::endl;
                std::cout << "Max iterations: " << M_maxIter << std::endl;
                std::cout << "Tolerance: " << M_tolerance << std::endl;
            }
            solveIterative_GMRES(M_precondType, M_maxIter, M_tolerance, M_verbose);
            break;
            
        case SolverType::BICGSTAB:
            if (M_verbose)
            {
                std::cout << "Solver: BiCGSTAB (iterative)" << std::endl;
                std::cout << "Max iterations: " << M_maxIter << std::endl;
                std::cout << "Tolerance: " << M_tolerance << std::endl;
            }
            solveIterative_BiCGSTAB(M_precondType, M_maxIter, M_tolerance, M_verbose);
            break;
            
        default:
            std::cerr << "ERROR: Unknown solver type, using SuperLU" << std::endl;
            solveDirect_SuperLU();
    }
    
    if (M_verbose)
    {
        if (M_converged)
        {
            std::cout << "✓ Solve successful";
            if (M_iterations > 0)
                std::cout << " (" << M_iterations << " iterations, residual: " 
                          << M_residual << ")";
            std::cout << std::endl;
        }
        else
        {
            std::cout << "✗ Solve failed or did not converge" << std::endl;
        }
        std::cout << "==========================\n" << std::endl;
    }
}


void LinearSystem::solveDirect_SuperLU()
{
    if (M_gotInverse)
    {
        gmm::clear(*M_Sol);
        gmm::mult(*M_InverseMatrix, *M_RHS, *M_Sol);
        M_converged = true;
    }
    else
    {
        scalar_type rcond;
        gmm::clear(*M_Sol);
        
        try
        {
            SuperLU_solve(*M_Matrix, *M_Sol, *M_RHS, rcond);
            M_converged = true;
            
            if (rcond < 1e-12)
            {
                std::cout << "WARNING: Matrix is near-singular (rcond = " 
                          << rcond << ")" << std::endl;
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "ERROR in SuperLU_solve: " << e.what() << std::endl;
            M_converged = false;
        }
    }
}

void LinearSystem::solveIterative_GMRES(PreconditionerType precondType,
                                       int maxIter,
                                       scalar_type tol,
                                       bool verbose)
{
    gmm::clear(*M_Sol);
    
    gmm::iteration iter(tol);
    iter.set_maxiter(maxIter);
    iter.set_noisy(verbose ? 1 : 0);
    
    try
    {
        switch (precondType)
        {
            case PreconditionerType::NONE:
            {
                if (verbose) std::cout << "Preconditioner: None" << std::endl;
                gmm::identity_matrix P;
                gmm::gmres(*M_Matrix, *M_Sol, *M_RHS, P, 50, iter);
                break;
            }
            
            case PreconditionerType::DIAGONAL:
            {
                if (verbose) std::cout << "Preconditioner: Diagonal (Jacobi)" << std::endl;
                gmm::diagonal_precond<sparseMatrix_Type> P(*M_Matrix);
                gmm::gmres(*M_Matrix, *M_Sol, *M_RHS, P, 50, iter);
                break;
            }
            
            case PreconditionerType::ILU:
            {
                if (verbose) std::cout << "Preconditioner: ILU(0)" << std::endl;
                gmm::ilu_precond<sparseMatrix_Type> P(*M_Matrix);
                gmm::gmres(*M_Matrix, *M_Sol, *M_RHS, P, 50, iter);
                break;
            }
            
            case PreconditionerType::ILUT:
            {
                if (verbose) std::cout << "Preconditioner: ILUT" << std::endl;
                gmm::ilut_precond<sparseMatrix_Type> P(*M_Matrix, 20, 1e-6);
                gmm::gmres(*M_Matrix, *M_Sol, *M_RHS, P, 50, iter);
                break;
            }
            
            default:
                std::cerr << "WARNING: Unknown preconditioner, using none" << std::endl;
                gmm::identity_matrix P;
                gmm::gmres(*M_Matrix, *M_Sol, *M_RHS, P, 50, iter);
        }
        
        M_iterations = iter.get_iteration();
        M_residual = iter.get_res();
        M_converged = iter.converged();
        
        if (!M_converged)
        {
            std::cerr << "WARNING: GMRES did not converge in " << maxIter 
                      << " iterations (residual: " << M_residual << ")" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "ERROR in GMRES solve: " << e.what() << std::endl;
        M_converged = false;
    }
}

void LinearSystem::solveIterative_BiCGSTAB(PreconditionerType precondType,
                                           int maxIter,
                                           scalar_type tol,
                                           bool verbose)
{
    gmm::clear(*M_Sol);
    
    gmm::iteration iter(tol);
    iter.set_maxiter(maxIter);
    iter.set_noisy(verbose ? 1 : 0);
    
    try
    {
        switch (precondType)
        {
            case PreconditionerType::NONE:
            {
                if (verbose) std::cout << "Preconditioner: None" << std::endl;
                gmm::identity_matrix P;
                gmm::bicgstab(*M_Matrix, *M_Sol, *M_RHS, P, iter);
                break;
            }
            
            case PreconditionerType::DIAGONAL:
            {
                if (verbose) std::cout << "Preconditioner: Diagonal (Jacobi)" << std::endl;
                gmm::diagonal_precond<sparseMatrix_Type> P(*M_Matrix);
                gmm::bicgstab(*M_Matrix, *M_Sol, *M_RHS, P, iter);
                break;
            }
            
            case PreconditionerType::ILU:
            {
                if (verbose) std::cout << "Preconditioner: ILU(0)" << std::endl;
                gmm::ilu_precond<sparseMatrix_Type> P(*M_Matrix);
                gmm::bicgstab(*M_Matrix, *M_Sol, *M_RHS, P, iter);
                break;
            }
            
            case PreconditionerType::ILUT:
            {
                if (verbose) std::cout << "Preconditioner: ILUT" << std::endl;
                gmm::ilut_precond<sparseMatrix_Type> P(*M_Matrix, 20, 1e-6);
                gmm::bicgstab(*M_Matrix, *M_Sol, *M_RHS, P, iter);
                break;
            }
            
            default:
                std::cerr << "WARNING: Unknown preconditioner, using none" << std::endl;
                gmm::identity_matrix P;
                gmm::bicgstab(*M_Matrix, *M_Sol, *M_RHS, P, iter);
        }
        
        M_iterations = iter.get_iteration();
        M_residual = iter.get_res();
        M_converged = iter.converged();
        
        if (!M_converged)
        {
            std::cerr << "WARNING: BiCGSTAB did not converge in " << maxIter 
                      << " iterations (residual: " << M_residual << ")" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "ERROR in BiCGSTAB solve: " << e.what() << std::endl;
        M_converged = false;
    }
}

#ifdef USE_MUMPS
void LinearSystem::solveDirect_MUMPS()
{
    if (M_verbose)
        std::cout << "[MUMPS] Starting parallel solve..." << std::endl;
    
    gmm::clear(*M_Sol);
    
    try {
        // MUMPS automatically uses MPI for parallel factorization
        gmm::MUMPS_solve(*M_Matrix, *M_Sol, *M_RHS);
        M_converged = true;
        
        if (M_verbose)
            std::cout << "[MUMPS] Solve completed successfully" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "[MUMPS] Error: " << e.what() << std::endl;
        M_converged = false;
        
        // Fallback to SuperLU
        std::cerr << "[MUMPS] Falling back to SuperLU..." << std::endl;
        solveDirect_SuperLU();
    }
}
#endif

void LinearSystem::computeInverse()
{
    std::vector<scalar_type> RHS(M_ndof,0.0);
    std::vector<scalar_type> provv(M_ndof,0.0);
    M_InverseMatrix.reset(new sparseMatrix_Type (M_ndof,M_ndof));

    for (size_type i=0; i<M_ndof;++i)
    {
	gmm::clear(RHS);
	gmm::clear(provv);

	RHS[i]=1;
        scalar_type rcond;
        SuperLU_solve(*M_Matrix, provv, RHS, rcond);
	gmm::copy(gmm::col_vector(provv), gmm::sub_matrix(*M_InverseMatrix, gmm::sub_interval(0, M_ndof),  gmm::sub_interval(0,1)));
	
    }
  
    M_gotInverse=true;

}

void LinearSystem::saveMatrix(const char* nomefile)
{
	gmm::MatrixMarket_IO::write(nomefile , *M_Matrix); 
}