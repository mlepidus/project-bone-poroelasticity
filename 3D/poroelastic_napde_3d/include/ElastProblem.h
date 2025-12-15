// ElastProblem.h - Quasi-static elasticity problem
#ifndef ELASTPROBLEM_H
#define ELASTPROBLEM_H

#include "Core.h"
#include "LinearSystem.h"
#include "ElastOperators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "TimeLoop.h"

/**
 * @class ElastProblem
 * @brief Quasi-static linear elasticity problem solver
 * 
 * Solves displacement field under applied loads and boundary conditions.
 * Supports poroelastic coupling through fluid pressure coupling term.
 */
class ElastProblem {
public:
    /**
     * @brief Construct elasticity problem from data file
     * @param dataFile GetPot parameter file
     * @param bulk Pointer to bulk domain
     */
    ElastProblem(const GetPot& dataFile, Bulk* bulk = NULL);
    
    /**
     * @brief Set time pointer for time-dependent BC and loads
     * @param timePtr Pointer to time loop manager
     */
    inline void setTimePointer(TimeLoop* timePtr) {
        M_time = timePtr;
    }
    
    /**
     * @brief Get finite element space for specified variable
     * @param where Domain selector
     * @param variable Variable name (e.g., "Displacement", "Stress")
     * @return Pointer to FEM object
     */
    FEM* getFEM(const std::string where = "bulk", const std::string variable = "Pressure");
    
    /// Get pointer to boundary condition handler
    inline BC* getBC() { return &M_BC; }
    
    /// Register DOFs with global linear system
    void addToSys(LinearSystem* sys);
    
    /**
     * @brief Assemble system matrix
     * @param sys Linear system object
     * @param where Domain selector
     */
    void assembleMatrix(LinearSystem* sys, std::string where);
    
    /**
     * @brief Assemble right-hand side vector
     * @param sys Linear system object
     * @param where Domain selector
     */
    void assembleRHS(LinearSystem* sys, std::string where);
    
    /**
     * @brief Enforce essential boundary conditions
     * @param firstTime If true, modify matrix; if false, only modify RHS
     */
    void enforceStrongBC(bool firstTime);
    
    /**
     * @brief Compute error against exact solution
     * @param what Variable to check
     * @param time Current time
     * @return L2 error norm
     */
    scalar_type computeError(std::string what, scalar_type time);
    
    /// Get pointer to linear system
    inline LinearSystem* getSys() { return M_Sys; }
    
    /// Solve linear system
    void solve();
    
    /// Update solution (old <- current)
    void updateSol();
    
    /// Extract solution from global system vector
    void extractSol(scalarVectorPtr_Type sol);
    
    /**
     * @brief Export solution to VTK format
     * @param folder Output directory
     * @param what Variables to export
     * @param frame Time frame number
     */
    void exportVtk(std::string folder = "./vtk", std::string what = "all", int frame = -1);
    
    /**
     * @brief Get number of degrees of freedom
     * @param variable Variable selector
     * @return Number of DOFs
     */
    size_type getNDOF(std::string variable = "all");
    
    /// Get previous time step solution
    scalarVectorPtr_Type getOldSol() { return M_DispSolOld; }
    
    /// Initialize problem (setup FEM spaces, matrices)
    void initialize();
    
    /// Compute normal stress on fault surface
    scalarVector_Type getNormalStressOnFault();
    
    /// Get friction coefficient
    inline scalar_type frictionCoeff() { return M_staticMu; }

private:
    TimeLoop* M_time;  ///< Pointer to time manager
    Bulk* M_Bulk;      ///< Pointer to bulk domain
    BC M_BC;           ///< Boundary condition handler
    
    // Finite element spaces
    FEM M_DispFEM;           ///< Displacement (vector) FEM
    FEM M_DispScalarFEM;     ///< Scalar displacement component FEM
    FEM M_CoeffFEM;          ///< Coefficient FEM for material properties
    FEM M_StressFEM;         ///< Stress (tensor) FEM
    FEM M_StressScalarFEM;   ///< Scalar stress component FEM
    
    LinearSystem* M_Sys;  ///< Pointer to global linear system
    
    // Solution vectors
    scalarVectorPtr_Type M_DispSol;                   ///< Current displacement
    scalarVectorPtr_Type M_DispSolOld;                ///< Previous displacement
    scalarVectorPtr_Type M_slipVel;                   ///< Slip velocity (fault)
    scalarVectorPtr_Type M_NormalStressSol;           ///< Normal stress vector
    scalarVectorPtr_Type M_NormalStressScalarSol;     ///< Normal stress scalar
    scalarVectorPtr_Type M_TangentStressScalarSol;    ///< Tangent stress scalar
    scalarVectorPtr_Type M_ST;                        ///< Tangential stress
    scalarVectorPtr_Type M_Slip;                      ///< Current slip
    scalarVectorPtr_Type M_SlipOld;                   ///< Previous slip
    
    // Auxiliary matrices
    sparseMatrixPtr_Type M_normalStressMat;
    sparseMatrixPtr_Type M_normalStressScalarMat;
    sparseMatrixPtr_Type M_tangentStressScalarMat;
    sparseMatrixPtr_Type M_massMatrixVector;
    sparseMatrixPtr_Type M_massMatrix;
    sparseMatrixPtr_Type M_massP2onFault;
    sparseMatrixPtr_Type M_maskUzawa;
    sparseMatrixPtr_Type M_glueMatrix;
    sparseMatrixPtr_Type M_normalGlueMatrix;
    sparseMatrixPtr_Type M_P02P2interp;
    sparseMatrixPtr_Type M_dispMassMatrix;
    
    size_type M_nbTotDOF;      ///< Total DOFs
    size_type M_nbTotBulkDOF;  ///< Bulk domain DOFs
    
    // Boundary condition data
    std::vector<size_type> M_rowsStrongBC;       ///< DOF indices with strong BC
    std::vector<size_type> M_rowsStrongBCFlags;  ///< Flags for BC type
    
    scalar_type M_staticMu;  ///< Static friction coefficient
    
    bool M_precomputeInterpolation;  ///< Flag for interpolation precomputation
    
    getfem::mesh_im M_intMethod;  ///< Integration method
    
    bool M_IsNitSym;   ///< Nitsche symmetry flag
    bool M_IsNitCons;  ///< Nitsche consistency flag
};

#endif // ELASTPROBLEM_H