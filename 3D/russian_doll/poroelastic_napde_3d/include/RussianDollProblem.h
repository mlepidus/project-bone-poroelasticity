// ============================================================================
// RussianDollProblem.h - Dual-porosity (PV + PLC) coupled problem
// ============================================================================
#ifndef RUSSIANDOLLPROBLEM_H
#define RUSSIANDOLLPROBLEM_H

#include "Core.h"
#include "LinearSystem.h"
#include "CoupledProblem.h"
//#include "PLCProblem.h"
#include "InterpolationManager.h"
#include "TimeLoop.h"

/**
 * @class RussianDollProblem
 * @brief Manages dual-porosity coupling between PV (vascular) and PLC (lacuno-canalicular)
 * 
 * The coupling occurs through leakage terms:
 * - PV equation gets: +γM p_l^N  (interpolated from PLC)
 * - PLC equation gets: +γM p_v^N (interpolated from PV)
 * 
 * Two coupling approaches are supported:
 * 1. LINE_INTERPOLATION: Extract PV solution along 1D lines, fit polynomial, apply as BC
 * 2. MESH_INTERPOLATION: Direct 3D mesh-to-mesh interpolation via GetFEM
 */

enum class CouplingApproach {
    LINE_INTERPOLATION,    // 1D line extraction + polynomial fitting
    MESH_INTERPOLATION     // Direct 3D-to-3D interpolation
};

class RussianDollProblem {
public:
    /**
     * @brief Construct from data file
     * @param dataFile GetPot parameter file
     * @param bulkPV Pointer to PV (outer) bulk domain
     * @param bulkPLC Pointer to PLC (inner) bulk domain  
     * @param time Pointer to time loop manager
     */
    RussianDollProblem(const GetPot& dataFile, 
                       Bulk* bulkPV, 
                       Bulk* bulkPLC,
                       TimeLoop* time);
    
    /// Destructor
    ~RussianDollProblem() = default;
    
    // ========================================================================
    // Setup Methods
    // ========================================================================
    
    /**
     * @brief Register the PV (vascular porosity) coupled problem
     * @param pvProblem Pointer to the PV CoupledProblem
     */
    void setPVProblem(CoupledProblem* pvProblem);
    
    /**
     * @brief Register the PLC (lacuno-canalicular) problem
     * @param plcProblem Pointer to the PLCProblem
     */
    void setPLCProblem(PLCProblem* plcProblem);
    
    /**
     * @brief Initialize both problems and interpolation structures
     */
    void initialize();
    
    // ========================================================================
    // Solution Methods - Staggered Coupling
    // ========================================================================
    
    /**
     * @brief Perform one time step with staggered coupling
     * 
     * Algorithm:
     * 1. Solve PV problem with p_l from previous step
     * 2. Interpolate p_v to PLC mesh
     * 3. Solve PLC problem with interpolated p_v
     * 4. (Optional) Iterate for convergence
     */
    void solveTimeStep();
    
    /**
     * @brief Solve PV problem for current time step
     */
    void solvePV();
    
    /**
     * @brief Solve PLC problem for current time step
     */
    void solvePLC();
    
    /**
     * @brief Transfer solution from PV to PLC domain
     * Uses the selected coupling approach
     */
    void interpolatePVtoPLC();
    
    /**
     * @brief Transfer solution from PLC to PV domain (for feedback)
     * Only needed if PLC feeds back to PV
     */
    void interpolatePLCtoPV();
    
    /**
     * @brief Update solutions for next time step
     */
    void updateSolutions();
    
    // ========================================================================
    // Assembly Methods
    // ========================================================================
    
    /**
     * @brief Assemble leakage coupling term for PV equation
     * Adds: +γ * ∫ φ_i * p_l dΩ to PV RHS
     * @param p_l_on_PV Pressure from PLC interpolated onto PV mesh
     */
    void assembleLeakageRHS_PV(const scalarVectorPtr_Type& p_l_on_PV);
    
    /**
     * @brief Assemble leakage coupling term for PLC equation  
     * Adds: +γ * ∫ φ_i * p_v dΩ to PLC RHS
     * @param p_v_on_PLC Pressure from PV interpolated onto PLC mesh
     */
    void assembleLeakageRHS_PLC(const scalarVectorPtr_Type& p_v_on_PLC);
    
    // ========================================================================
    // Output Methods
    // ========================================================================
    
    /**
     * @brief Export solutions to VTK
     * @param folder Output directory
     * @param frame Time frame number
     */
    void exportVtk(const std::string& folder, int frame);
    
    /**
     * @brief Compute errors against exact solutions (if available)
     * @param time Current time
     * @return Vector of errors [p_v_error, u_v_error, p_l_error, u_l_error]
     */
    std::vector<scalar_type> computeErrors(scalar_type time);
    
    // ========================================================================
    // Getters
    // ========================================================================
    
    inline CoupledProblem* getPVProblem() { return M_pvProblem; }
    inline PLCProblem* getPLCProblem() { return M_plcProblem; }
    inline InterpolationManager* getInterpolationManager() { return M_interpManager.get(); }
    
    /// Get leakage coefficient γ
    inline scalar_type getLeakageCoeff() const { return M_gamma; }

private:
    // Domain pointers
    Bulk* M_bulkPV;           ///< PV (outer) bulk domain
    Bulk* M_bulkPLC;          ///< PLC (inner) bulk domain
    TimeLoop* M_time;         ///< Time manager
    
    // Sub-problems
    CoupledProblem* M_pvProblem;    ///< PV poroelasticity problem
    //PLCProblem* M_plcProblem;       ///< PLC poroelasticity problem
    
    // Interpolation
    std::unique_ptr<InterpolationManager> M_interpManager;
    CouplingApproach M_couplingApproach;
    
    // Coupling parameters
    scalar_type M_gamma;             ///< Leakage coefficient
    size_type M_maxCouplingIter;     ///< Max iterations for staggered coupling
    scalar_type M_couplingTolerance; ///< Convergence tolerance
    
    // Solution transfer vectors
    scalarVectorPtr_Type M_pv_on_PLC;  ///< PV pressure interpolated to PLC mesh
    scalarVectorPtr_Type M_pl_on_PV;   ///< PLC pressure interpolated to PV mesh
    
    // Previous step solutions (for coupling)
    scalarVectorPtr_Type M_pv_old;     ///< Previous PV pressure
    scalarVectorPtr_Type M_pl_old;     ///< Previous PLC pressure
    
    // Configuration
    std::string M_section;           ///< Data file section
    bool M_initialized;              ///< Initialization flag
};

#endif // RUSSIANDOLLPROBLEM_H