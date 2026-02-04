// CoupledProblem.h - Strongly coupled poroelasticity problem
#ifndef COUPLEDPROBLEM_H
#define COUPLEDPROBLEM_H

#include "Core.h"
#include "LinearSystem.h"
#include "ElastProblem.h"
#include "DarcyProblemT.h"
#include "StringUtility.h"
#include "TimeLoop.h"

/**
 * @class CoupledProblem
 * @brief Fully coupled poroelasticity problem (Biot's equations)
 * 
 * Manages strong coupling between Darcy flow and linear elasticity.
 * Assembles monolithic system with pressure-stress coupling terms.
 * Note: Currently supports linear slip only.
 */
class CoupledProblem {
public:
    /**
     * @brief Construct coupled problem from data file
     * @param dataFile GetPot parameter file
     * @param bulk Pointer to bulk domain
     * @param time Pointer to time loop manager
     * @param section Data file section (default: root, can be "plc/" for PLCProblem)
     */
    CoupledProblem(const GetPot& dataFile, Bulk* bulk = NULL, 
                   TimeLoop* time = NULL, const std::string& section = "");
    
    /**
     * @brief Virtual destructor for proper inheritance
     */
    virtual ~CoupledProblem() = default;

    /**
     * @brief Register elasticity sub-problem
     * @param elast Pointer to ElastProblem object
     */
    inline void addElastPB(ElastProblem* elast) {
        M_ElastPB = elast;
        M_ElastPB->setTimePointer(M_time);
    }
    
    /**
     * @brief Register Darcy flow sub-problem
     * @param darcy Pointer to DarcyProblemT object
     */
    inline void addDarcyPB(DarcyProblemT* darcy) {
        M_DarcyPB = darcy;
        M_DarcyPB->setTimePointer(M_time);
    }
    
    /**
     * @brief Compute solution error against exact solution
     * @param t Current time
     * @return Error vector (pressure error, displacement error)
     */
    bgeot::base_node computeError(scalar_type t); 
    
    
    /// Register DOFs with global linear system
    void addToSys(LinearSystem* sys); 
    
    /**
     * @brief Assemble global system matrix
     */
    virtual void assembleMatrix(); 
    
    /**
     * @brief Assemble global right-hand side vector
     */
    virtual void assembleRHS();
    
    /// Add sub-problem matrices to global system
    void addSubSystems();
    
    /// Add sub-problem RHS vectors to global system
    void addSubSystemsRHS();
    
    /// Clear sub-system matrices
    void clearSubSystems();
    
    /// Clear sub-system RHS vectors
    void clearSubSystemsRHS();
    
    /// Solve coupled system
    virtual void solve();
    
    /// Update solution for next time step
    virtual void updateSol();
    
    /**
     * @brief Export solution to VTK format
     * @param folder Output directory
     * @param what Variables to export ("all", "darcy", "elast")
     * @param frame Time frame number
     */
    virtual void exportVtk(std::string folder = "./vtk", std::string what = "all", int frame = -1);
    
    /**
     * @brief Export time history data to file
     * @param timestepData Data for current time step
     * @param filename Output file name
     * @param step Time step number
     */
    void exportHistory(const scalarVector_Type& timestepData,
                      const std::string& filename, size_t step);
    
    /**
     * @brief Enforce boundary conditions
     * @param firstTime If true, modify matrix and RHS; if false, only RHS
     */
    virtual void enforceStrongBC(bool firstTime);

    // ========================================================================
    // Getters (for PLCProblem and RussianDollProblem access)
    // ========================================================================
    
    /// Get pointer to elasticity problem
    inline ElastProblem* getElastPB() { return M_ElastPB; }
    
    /// Get pointer to Darcy problem
    inline DarcyProblemT* getDarcyPB() { return M_DarcyPB; }
    
    /// Get pointer to linear system
    inline LinearSystem* getSys() { return M_Sys; }
    
    /// Get pointer to bulk domain
    inline Bulk* getBulk() { return M_Bulk; }
    
    /// Get pointer to time manager
    inline TimeLoop* getTime() { return M_time; }
    
    /// Get pressure solution from Darcy sub-problem
    inline scalarVectorPtr_Type getPressure() {
        if (M_DarcyPB) {
            scalarVectorPtr_Type pressure;
            pressure.reset(new scalarVector_Type(M_DarcyPB->getPressureSolution()));
            return pressure;
        }
        return scalarVectorPtr_Type();
    }
    
    /// Get displacement solution
    inline scalarVectorPtr_Type getDisplacement() { return M_solEl; }
    
    /// Get pressure mesh_fem
    inline const getfem::mesh_fem& getMfPressure() const {
        return *(M_DarcyPB->getFEM("Pressure")->getFEM());
    }
    
    /// Get displacement mesh_fem
    inline const getfem::mesh_fem& getMfDisplacement() const {
        return *(M_ElastPB->getFEM()->getFEM());
    }
    
    /// Get total DOFs
    inline size_type getNbTotDOF() const { return M_nbTotDOF; }
    
    /// Get number of pressure DOFs
    size_type getNbPressureDOF() const { 
        return M_DarcyPB ? M_DarcyPB->getNDOF("Pressure") : 0; 
    }
    
    /// Get number of velocity DOFs
    size_type getNbVelocityDOF() const { 
        return M_DarcyPB ? M_DarcyPB->getNDOF("Velocity") : 0; 
    }
    
    /// Get number of elasticity DOFs
    size_type getNbElastDOF() const { 
        return M_ElastPB ? M_ElastPB->getNDOF() : 0; 
    }

protected:

    Bulk* M_Bulk;          ///< Pointer to bulk domain
    TimeLoop* M_time;      ///< Pointer to time manager
    std::string M_section; ///< Data file section prefix

    LinearSystem* M_Sys;       ///< Global coupled system
    LinearSystem M_elastSys;   ///< Elasticity sub-system
    LinearSystem M_darcySys;   ///< Darcy sub-system
    
    ElastProblem* M_ElastPB;   ///< Pointer to elasticity problem
    DarcyProblemT* M_DarcyPB;  ///< Pointer to Darcy problem
    
    scalarVectorPtr_Type M_solEl;  ///< Elasticity solution
    scalarVectorPtr_Type M_solDa;  ///< Darcy solution
    
    sparseMatrixPtr_Type M_PressureStress;  ///< Poroelastic coupling matrix
    
    getfem::mesh_im M_intMethod;  ///< Integration method
    
    size_type M_nbTotDOF;  ///< Total coupled system DOFs
    size_type step;        ///< Time step counter
    
    mutable LifeV::Parser M_parser;  ///< Expression parser
};

#endif // COUPLEDPROBLEM_H
