// ============================================================================
// LinearSystem.h - Linear system assembly and solution
// ============================================================================
#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "Core.h"

/**
 * @enum SolverType
 * @brief Available solver types for the linear system
 */
enum class SolverType {
    SUPERLU,      ///< Direct solver (SuperLU) - always works, memory intensive
    GMRES,        ///< Iterative GMRES - good for non-symmetric systems
    BICGSTAB      ///< Iterative BiCGSTAB - often faster than GMRES
};

/**
 * @enum PreconditionerType
 * @brief Available preconditioners for iterative solvers
 */
enum class PreconditionerType {
    NONE,         ///< No preconditioning
    DIAGONAL,     ///< Diagonal (Jacobi) preconditioner
    ILU,          ///< Incomplete LU factorization
    ILUT          ///< ILU with threshold
};


/**
 * @class LinearSystem
 * @brief Container for sparse linear system Ax = b with block assembly support
 * 
 * Manages system matrix, right-hand side vector, and solution vector.
 * Supports block-wise assembly for coupled multi-physics problems.
 * Solver configuration is read from data file and stored internally.
 */
class LinearSystem {
public:
    /// Default constructor
    LinearSystem();
    
    /**
     * @brief Initialize system matrix and vectors
     * @param ndof Total number of degrees of freedom
     */
    void addToMatrix(int ndof);
    
    /**
     * @brief Copy submatrix into system matrix (overwrites existing entries)
     * @param M Source submatrix
     * @param first_row Starting row index
     * @param first_column Starting column index
     * @param scale Scaling factor (default: 1.0)
     * @param transpose If true, transpose M during copy
     */
    void copySubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column,
                      scalar_type scale = 1.0, bool transpose = false);
    
    /**
     * @brief Add submatrix to system matrix (accumulates with existing entries)
     * @param M Source submatrix
     * @param first_row Starting row index
     * @param first_column Starting column index
     * @param scale Scaling factor (default: 1.0)
     * @param transpose If true, transpose M during addition
     */
    void addSubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column,
                     scalar_type scale = 1.0, bool transpose = false);
    
    /// Copy subvector into RHS (overwrites existing entries)
    void copySubVector(scalarVectorPtr_Type V, int first_row, scalar_type scale = 1.0);
    
    /// Add subvector to RHS (accumulates with existing entries)
    void addSubVector(scalarVectorPtr_Type V, int first_row, scalar_type scale = 1.0);
    
    /**
     * @brief Extract portion of solution or RHS vector
     * @param V Destination vector
     * @param first_row Starting index
     * @param where Source selector: "sol" for solution, "rhs" for right-hand side
     */
    void extractSubVector(scalarVectorPtr_Type V, int first_row, std::string where = "sol");
    
    /// Get pointer to system matrix
    inline sparseMatrixPtr_Type getMatrix() { return M_Matrix; }
    
    /// Get pointer to right-hand side vector
    inline scalarVectorPtr_Type getRHS() { return M_RHS; }
    
    /// Get pointer to solution vector
    inline scalarVectorPtr_Type getSol() { return M_Sol; }
    
    /**
     * @brief Add subsystem into block structure (matrix and RHS)
     * @param small Subsystem to add
     * @param shiftRows Row offset in global system
     * @param shiftColumns Column offset in global system
     */
    void addSubSystem(LinearSystem* small, size_type shiftRows, size_type shiftColumns);
    
    /**
     * @brief Add subsystem RHS only (for time-stepping with constant matrix)
     * @param small Subsystem to add
     * @param shiftRows Row offset in global system
     */
    void addSubSystemRHS(LinearSystem* small, size_type shiftRows);
    
    /**
     * @brief Configure solver from data file and store settings internally
     * @param dataFile GetPot parameter file
     * @param section Section name in data file (default: "solver/")
     * 
     * Reads solver configuration from data file and stores in member variables.
     * Expected format in data file:
     * [solver]
     *     type = BICGSTAB
     *     preconditioner = ILU
     *     maxIterations = 1000
     *     tolerance = 1.0e-8
     *     verbose = 1
     * [../]
     */
    void configureSolver(const GetPot& dataFile, const std::string& section = "solver/");
    
    /**
     * @brief Solve the linear system using stored configuration
     * Uses the solver settings loaded via configureSolver()
     */
    void solve();
    
    /// Compute and store matrix inverse (for small systems only)
    void computeInverse();
    
    /// Export system matrix to Matrix Market format file
    void saveMatrix(const char* filename = "Matrix.mm");
    
    /**
     * @brief Multiply submatrix block by vector and add to RHS
     * @param V Vector to multiply
     * @param first_row Starting row in system matrix
     * @param first_column Starting column in system matrix
     * @param nrows Number of rows to multiply
     * @param ncols Number of columns to multiply
     */
    void multAddToRHS(scalarVectorPtr_Type V, int first_row, int first_column,
                     size_type nrows, size_type ncols);
    
    /// Multiply arbitrary matrix M by vector V and add result to RHS
    void multAddToRHS(sparseMatrixPtr_Type M, scalarVectorPtr_Type V,
                     int first_rowVector, int first_rowRHS,
                     scalar_type scale = 1.0, bool transposed = false);
    
    /// Overload accepting non-pointer vector
    void multAddToRHS(sparseMatrixPtr_Type M, scalarVector_Type& V,
                     int first_rowVector, int first_rowRHS,
                     scalar_type scale = 1.0, bool transposed = false);
    
    /// Clear right-hand side vector (set all entries to zero)
    inline void cleanRHS() { gmm::clear(*M_RHS); }
    
    /// Clear system matrix (set all entries to zero)
    inline void cleanMAT() { gmm::clear(*M_Matrix); }
    
    /// Set entire row to zero (for applying constraints)
    inline void setNullRow(size_type which) {
        for (size_type j = 0; j < M_ndof; ++j)
            (*M_Matrix)(which, j) = 0;
    }
    
    /// Set individual matrix entry
    inline void setMatrixValue(size_type i, size_type j, scalar_type value) {
        (*M_Matrix)(i, j) = value;
    }
    
    /// Set individual RHS entry
    inline void setRHSValue(size_type i, scalar_type value) {
        (*M_RHS)[i] = value;
    }

    /// Get solver convergence information
    inline int getIterations() const { return M_iterations; }
    inline scalar_type getResidual() const { return M_residual; }
    inline bool hasConverged() const { return M_converged; }

private:
    sparseMatrixPtr_Type M_Matrix;        ///< System matrix A
    sparseMatrixPtr_Type M_InverseMatrix; ///< Inverse matrix A^{-1} (if computed)
    scalarVectorPtr_Type M_RHS;           ///< Right-hand side vector b
    scalarVectorPtr_Type M_Sol;           ///< Solution vector x

    bool M_gotInverse;                    ///< Flag indicating if inverse is available
    size_type M_ndof;                     ///< Number of degrees of freedom

    // Solver convergence information
    int M_iterations;                     ///< Number of iterations performed
    scalar_type M_residual;               ///< Final residual norm
    bool M_converged;                     ///< Convergence flag
    
    // Solver configuration (stored internally, read from data file)
    SolverType M_solverType;              ///< Selected solver type
    PreconditionerType M_precondType;     ///< Selected preconditioner type
    int M_maxIter;                        ///< Maximum iterations for iterative solvers
    scalar_type M_tolerance;              ///< Convergence tolerance
    bool M_verbose;                       ///< Verbose output flag
    
    // Private solver methods
    void solveDirect_SuperLU();
    void solveIterative_GMRES(PreconditionerType precond, int maxIter, 
                             scalar_type tol, bool verbose);
    void solveIterative_BiCGSTAB(PreconditionerType precond, int maxIter,
                                scalar_type tol, bool verbose);
};

#endif // LINEARSYSTEM_H