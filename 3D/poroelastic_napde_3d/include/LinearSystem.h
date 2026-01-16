// ============================================================================
// LinearSystem.h - Linear system assembly and solution
// ============================================================================
#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "Core.h"

enum class SolverType {
    SUPERLU,
    GMRES,
    BICGSTAB
};

enum class PreconditionerType {
    NONE,
    DIAGONAL,
    ILU,
    ILUT
};

class LinearSystem {
public:
    /// Default constructor (no solver configuration)
    LinearSystem();
    
    /// Constructor with automatic solver configuration from data file
    LinearSystem(const GetPot& dataFile, const std::string& section = "solver/");
    
    void addToMatrix(int ndof);
    
    void copySubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column,
                      scalar_type scale = 1.0, bool transpose = false);
    
    void addSubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column,
                     scalar_type scale = 1.0, bool transpose = false);
    
    void copySubVector(scalarVectorPtr_Type V, int first_row, scalar_type scale = 1.0);
    void addSubVector(scalarVectorPtr_Type V, int first_row, scalar_type scale = 1.0);
    void extractSubVector(scalarVectorPtr_Type V, int first_row, std::string where = "sol");
    
    inline sparseMatrixPtr_Type getMatrix() { return M_Matrix; }
    inline scalarVectorPtr_Type getRHS() { return M_RHS; }
    inline scalarVectorPtr_Type getSol() { return M_Sol; }
    
    void addSubSystem(LinearSystem* small, size_type shiftRows, size_type shiftColumns);
    void addSubSystemRHS(LinearSystem* small, size_type shiftRows);
    
    /**
     * @brief Configure solver from data file
     * Can be called manually if default constructor was used
     */
    void configureSolver(const GetPot& dataFile, const std::string& section = "solver/");
    
    /**
     * @brief Solve using stored configuration
     */
    void solve();
    
    // ========================================================================
    // NEW: Simple setter methods for runtime parameter changes
    // ========================================================================
    
    /// Change solver type at runtime
    void setSolverType(SolverType type) { 
        M_solverType = type; 
        if (M_verbose) std::cout << "Solver type changed" << std::endl;
    }
    
    /// Change preconditioner type
    void setPreconditioner(PreconditionerType type) { 
        M_precondType = type; 
        if (M_verbose) std::cout << "Preconditioner changed" << std::endl;
    }
    
    /// Set maximum iterations for iterative solvers
    void setMaxIterations(int maxIter) { 
        M_maxIter = maxIter; 
        if (M_verbose) std::cout << "Max iterations set to " << maxIter << std::endl;
    }
    
    /// Set convergence tolerance
    void setTolerance(scalar_type tol) { 
        M_tolerance = tol; 
        if (M_verbose) std::cout << "Tolerance set to " << tol << std::endl;
    }
    
    /// Enable/disable verbose output
    void setVerbose(bool verbose) { M_verbose = verbose; }
    
    /// Get current solver type
    SolverType getSolverType() const { return M_solverType; }
    
    /// Get convergence information
    inline int getIterations() const { return M_iterations; }
    inline scalar_type getResidual() const { return M_residual; }
    inline bool hasConverged() const { return M_converged; }
    
    void computeInverse();
    void saveMatrix(const char* filename = "Matrix.mm");
    
    void multAddToRHS(scalarVectorPtr_Type V, int first_row, int first_column,
                     size_type nrows, size_type ncols);
    
    void multAddToRHS(sparseMatrixPtr_Type M, scalarVectorPtr_Type V,
                     int first_rowVector, int first_rowRHS,
                     scalar_type scale = 1.0, bool transposed = false);
    
    void multAddToRHS(sparseMatrixPtr_Type M, scalarVector_Type& V,
                     int first_rowVector, int first_rowRHS,
                     scalar_type scale = 1.0, bool transposed = false);
    
    inline void cleanRHS() { gmm::clear(*M_RHS); }
    inline void cleanMAT() { gmm::clear(*M_Matrix); }
    
    inline void setNullRow(size_type which) {
        for (size_type j = 0; j < M_ndof; ++j)
            (*M_Matrix)(which, j) = 0;
    }
    
    inline void setMatrixValue(size_type i, size_type j, scalar_type value) {
        (*M_Matrix)(i, j) = value;
    }
    
    inline void setRHSValue(size_type i, scalar_type value) {
        (*M_RHS)[i] = value;
    }

private:
    sparseMatrixPtr_Type M_Matrix;
    sparseMatrixPtr_Type M_InverseMatrix;
    scalarVectorPtr_Type M_RHS;
    scalarVectorPtr_Type M_Sol;

    bool M_gotInverse;
    size_type M_ndof;

    // Convergence information
    int M_iterations;
    scalar_type M_residual;
    bool M_converged;
    
    // Solver configuration (stored internally)
    SolverType M_solverType;
    PreconditionerType M_precondType;
    int M_maxIter;
    scalar_type M_tolerance;
    bool M_verbose;
    bool M_configured;  // Flag to check if solver was configured
    
    // Private solver methods
    void solveDirect_SuperLU();
    void solveIterative_GMRES(PreconditionerType precond, int maxIter, 
                             scalar_type tol, bool verbose);
    void solveIterative_BiCGSTAB(PreconditionerType precond, int maxIter,
                                scalar_type tol, bool verbose);
};

#endif // LINEARSYSTEM_H