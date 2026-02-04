// ============================================================================
// Problem Class Headers
// ============================================================================

// DarcyProblemT.h - Time-dependent Darcy flow problem
#ifndef DARCYPROBLEMT_H
#define DARCYPROBLEMT_H

#include "Core.h"
#include "LinearSystem.h"
#include "DarcyOperators.h"
#include "StringUtility.h"
#include "TimeLoop.h"

#ifndef BOUNDARY_ASSIGNMENT_TYPE_DEFINED
#define BOUNDARY_ASSIGNMENT_TYPE_DEFINED
enum class BoundaryAssignmentType {
    GEOMETRIC_CYLINDER,      ///< Automatic detection based on face normals and position
    GEOMETRIC_SQUARE,        ///< Automatic detection for square geometry
    TAG_NAME,       ///< Match Gmsh physical names (outer, inner, top, bottom, etc.)
    TAG_NUMBER      ///< Select N largest/smallest tag numbers
};
#endif

/**
 * @class DarcyProblemT
 * @brief Time-dependent Darcy flow problem solver
 * 
 * Implements mixed finite element formulation for transient porous media
 * flow. Manages pressure and velocity fields with backward Euler time
 * discretization.
 */
class DarcyProblemT {
public:
    /**
     * @brief Construct Darcy problem from data file
     * @param dataFile GetPot parameter file
     * @param bulk Pointer to bulk domain
     */
    DarcyProblemT(const GetPot& dataFile, Bulk* bulk = NULL, const std::string& basePath = "bulkData/");
    
    /**
     * @brief Set time loop pointer for time-dependent BC and source terms
     * @param timePtr Pointer to time loop manager
     */
    inline void setTimePointer(TimeLoop* timePtr) {
        M_timeLoop = timePtr;
        M_dt = timePtr->dt();
    }
    
    /**
     * @brief Get finite element space for specified variable
     * @param where Domain selector ("bulk" or "fracture")
     * @param variable Variable name ("Pressure" or "Velocity")
     * @return Pointer to FEM object
     */
    FEM* getFEM(const std::string variable = "Pressure");
    
    /// Initialize problem (setup FEM spaces, initial conditions)
    void initialize();
    
    /// Register DOFs with global linear system
    void addToSys(LinearSystem* sys);
    
    /**
     * @brief Assemble system matrix
     * @param sys Linear system object
     * @param where Domain selector
     */
    void assembleMatrix();
    
    /**
     * @brief Assemble right-hand side vector
     * @param sys Linear system object
     * @param where Domain selector
     */
    void assembleRHS();
    
    /// Solve linear system for current time step
    void solve();
    
    /// Update solution (old <- current) for next time step
    void updateSol();
    
    /**
     * @brief Export solution to VTK format
     * @param folder Output directory
     * @param what Variables to export ("all", "pressure", "velocity")
     * @param frame Time frame number (-1 for steady state)
     */
    void exportVtk(std::string folder = "./vtk", std::string what = "all", int frame = -1);
    
    /// Get const reference to pressure solution vector
    const scalarVector_Type& getPressureSolution() const { return *M_pressureSol; }
    
    /// Get const reference to velocity solution vector
    const scalarVector_Type& getVelocitySolution() const { return *M_velocitySol; }
    
    /**
     * @brief Compute error against exact solution
     * @param what Variable to check ("all", "pressure", "velocity")
     * @param time Current time for exact solution evaluation
     * @return L2 error norm
     */
    scalar_type computeError(std::string what = "all", scalar_type time = 0);
    
    /**
     * @brief Get number of degrees of freedom
     * @param variable Variable selector ("all", "pressure", "velocity")
     * @return Number of DOFs
     */
    size_type getNDOF(std::string variable = "all");
    
    /// Extract solution from global system vector
    void extractSol(scalarVectorPtr_Type sol);

    /// Get BC object
    inline BC* getBC() { return &M_BC; }
    
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

    scalar_type M_dt;                       ///< Time step size
    TimeLoop* M_timeLoop;                   ///< Pointer to time manager
    Bulk* M_Bulk;                           ///< Pointer to bulk domain
    BC M_BC;                                ///< Boundary condition handler
    
    // Boundary assignment configuration
    BoundaryAssignmentType M_boundaryType;  ///< How to assign boundary regions
    bool M_tagNumberLargest;                ///< For TAG_NUMBER: select largest (true) or smallest (false)
    
    // Finite element spaces
    FEM M_PressureFEM;        ///< Pressure (scalar) FEM
    FEM M_CoeffFEM;           ///< Coefficient FEM for material properties
    FEM M_VelocityFEM;        ///< Velocity (H(div)) FEM
    FEM M_visualizationVFEM;  ///< Velocity FEM for visualization
    
    LinearSystem* M_Sys;  ///< Pointer to global linear system
    
    // Solution vectors
    scalarVectorPtr_Type M_pressureSol;     ///< Current pressure
    scalarVectorPtr_Type M_pressureSolOld;  ///< Previous time step pressure
    scalarVectorPtr_Type M_pressureSolIni;  ///< Initial pressure
    scalarVectorPtr_Type M_velocitySol;     ///< Current velocity
    
    sparseMatrixPtr_Type M_pressureMass;  ///< Pressure-pressure block matrix

    size_type M_nbTotDOF;      ///< Total DOFs (pressure + velocity)
    size_type M_nbTotBulkDOF;  ///< Bulk domain DOFs only
    
    getfem::mesh_im M_intMethod;  ///< Integration method
    
    mutable LifeV::Parser M_parser;  ///< Expression parser
};

#endif // DARCYPROBLEMT_H