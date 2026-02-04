// ElastProblem.h - Quasi-static linear elasticity with poroelastic coupling
#ifndef ELASTPROBLEM_H
#define ELASTPROBLEM_H

#include "Core.h"
#include "LinearSystem.h"
#include "ElastOperators.h"
#include "StringUtility.h"
#include "TimeLoop.h"

/**
 * @enum BoundaryAssignmentType
 * @brief Specifies how boundary regions are assigned from mesh data
 */
#ifndef BOUNDARY_ASSIGNMENT_TYPE_DEFINED
#define BOUNDARY_ASSIGNMENT_TYPE_DEFINED
enum class BoundaryAssignmentType {
    GEOMETRIC_CYLINDER,      ///< Automatic detection based on face normals and position (CENTERED ON (1,1))
    GEOMETRIC_SQUARE,        ///< Automatic detection for square geometry
    TAG_NAME,       ///< Match Gmsh physical names (outer, inner, top, bottom, etc.)
    TAG_NUMBER      ///< Select N largest/smallest tag numbers
};
#endif

/**
 * @class ElastProblem
 * @brief Quasi-static linear elasticity problem solver with poroelastic coupling
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
    ElastProblem(const GetPot& dataFile, Bulk* bul=NULL, 
                          const std::string& basePath = "bulkData/");
    
    /**
     * @brief Set time pointer for time-dependent BC and loads
     * @param timePtr Pointer to time loop manager
     */
    inline void setTimePointer(TimeLoop* timePtr) {
        M_time = timePtr;
    }
    
    /**
     * @brief Get finite element space for displacement
     * @return Pointer to FEM object
     */
    FEM* getFEM() { return &M_DispFEM; }
    
    /// Get pointer to boundary condition handler
    inline BC* getBC() { return &M_BC; }
    
    /// Register DOFs with global linear system
    void addToSys(LinearSystem* sys);
    
    /**
     * @brief Assemble system matrix
     */
    void assembleMatrix();
    
    /**
     * @brief Assemble right-hand side vector
     */
    void assembleRHS();
    
    /**
     * @brief Enforce essential boundary conditions
     * @param firstTime If true, modify matrix; if false, only modify RHS
     */
    void enforceStrongBC(bool firstTime);
    
    /**
     * @brief Compute error against exact solution
     * @param what Variable to check (only "Displacement" supported)
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
     * @param what Variables to export (only "Displacement" supported)
     * @param frame Time frame number
     */
    void exportVtk(std::string folder = "./vtk", int frame = -1);
    
    /**
     * @brief Get number of degrees of freedom
     * @param variable Variable selector (only "Disp" or "all" supported)
     * @return Number of DOFs
     */
    size_type getNDOF(std::string variable = "all");
    
    /// Get previous time step solution
    scalarVectorPtr_Type getOldSol() { return M_DispSolOld; }
    
    /// Initialize problem (setup FEM spaces, matrices)
    void initialize();

private:
    /**
     * @brief Setup boundary conditions based on assignment type
     */
    void setupBoundaryConditions();
    
    /**
     * @brief Parse boundary assignment type from string
     * @param typeStr String from input file ("geometric", "tagName", "tagNumber")
     * @return Corresponding enum value
     */
    static BoundaryAssignmentType parseBoundaryType(const std::string& typeStr);

    TimeLoop* M_time;      ///< Pointer to time manager
    Bulk* M_Bulk;          ///< Pointer to bulk domain
    BC M_BC;               ///< Boundary condition handler
    
    // Boundary assignment configuration
    BoundaryAssignmentType M_boundaryType;  ///< How to assign boundary regions
    bool M_tagNumberLargest;                ///< For TAG_NUMBER: select largest (true) or smallest (false)
    
    // Finite element spaces
    FEM M_DispFEM;         ///< Displacement (vector) FEM
    FEM M_CoeffFEM;        ///< Coefficient FEM for material properties
    
    LinearSystem* M_Sys;   ///< Pointer to global linear system
    
    // Solution vectors
    scalarVectorPtr_Type M_DispSol;     ///< Current displacement
    scalarVectorPtr_Type M_DispSolOld;  ///< Previous displacement
    
    // Mass matrix for error computation
    sparseMatrixPtr_Type M_dispMassMatrix;
    
    size_type M_nbTotDOF;  ///< Total DOFs
    
    // Boundary condition data for strong enforcement
    std::vector<size_type> M_rowsStrongBC;       ///< DOF indices with strong BC
    std::vector<size_type> M_rowsStrongBCFlags;  ///< Flags for BC type
    
    getfem::mesh_im M_intMethod;  ///< Integration method
};

#endif // ELASTPROBLEM_H