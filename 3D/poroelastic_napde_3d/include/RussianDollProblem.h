// ============================================================================
// RussianDollProblem.h - Dual-porosity (PV + PLC) coupled problem
// ============================================================================
// This version incorporates the InterpolationManager for one-way coupling
// from PV to PLC with polynomial BC coefficients.
// ============================================================================
#ifndef RUSSIANDOLLPROBLEM_H
#define RUSSIANDOLLPROBLEM_H

#include "Core.h"
#include "LinearSystem.h"
#include "CoupledProblem.h"
#include "PLCProblem.h"
#include "InterpolationManager.h"
#include "TimeLoop.h"

/**
 * @enum CouplingApproach
 * @brief Defines the method used to transfer pressure from PV to PLC
 */
enum class CouplingApproach {
    LINE_INTERPOLATION,    ///< 1D line extraction + polynomial fitting
    MESH_INTERPOLATION     ///< Direct 3D-to-3D interpolation via GetFEM
};

/**
 * @class RussianDollProblem
 * @brief Manages dual-porosity coupling between PV (vascular) and PLC (lacuno-canalicular)
 * 
 * This class orchestrates a one-way coupling from PV to PLC:
 * - PV is solved first (uncoupled from PLC, i.e., p_l = 0 in PV equation)
 * - PV pressure is interpolated to PLC boundary via line extraction + polynomial fit
 * - PLC is solved with the interpolated PV pressure as Dirichlet BC on outer wall
 * 
 * The governing equations are:
 * 
 * PV (Vascular Porosity):
 *   K u_v^N - α_v C p_v^N = R_v^N
 *   α_v C^T u̇_v^N + H_v p_v^N + (1/M_v) M ṗ_v^N = Q_v^N
 * 
 * PLC (Lacuno-canalicular Porosity):
 *   K u_l^N - α_l C p_l^N = R_l^N
 *   α_l C^T u̇_l^N + H_l p_l^N + (1/M_l) M ṗ_l^N = Q_l^N
 *   with BC: p_l = p_v(z) on outer wall (Dirichlet from interpolated PV pressure)
 * 
 * Key Features:
 * - InterpolationManager is incorporated directly for clean coefficient management
 * - Polynomial coefficients from PV → PLC transfer are stored for BC evaluation
 * - PLCProblem's BC class can access these coefficients through callbacks
 * 
 * Typical Usage:
 * @code
 *   RussianDollProblem russian(dataFile, bulkPV, bulkPLC, time);
 *   russian.setPVProblem(&pvProblem);
 *   russian.setPLCProblem(&plcProblem);
 *   russian.initialize();
 *   
 *   for each time step:
 *       russian.solveTimeStep();
 *       russian.exportVtk(folder, step);
 *       russian.updateSolutions();
 * @endcode
 */
class RussianDollProblem {
public:
    // ========================================================================
    // Construction / Destruction
    // ========================================================================
    
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
     * 
     * Also sets up the BC callback for polynomial-based Dirichlet conditions
     */
    void setPLCProblem(PLCProblem* plcProblem);
    
    /**
     * @brief Initialize both problems and interpolation structures
     * 
     * This must be called after setPVProblem() and setPLCProblem()
     */
    void initialize();
    
    // ========================================================================
    // Solution Methods - Staggered One-Way Coupling
    // ========================================================================
    
    /**
     * @brief Perform one time step with one-way coupling (PV → PLC)
     * 
     * Algorithm:
     * 1. Solve PV problem (standalone, p_l = 0)
     * 2. Extract PV pressure along vertical line
     * 3. Fit polynomial p_v(z)
     * 4. Update PLC boundary coefficients
     * 5. Solve PLC problem with polynomial BC on outer wall
     */
    void solveTimeStep();
    
    /**
     * @brief Solve PV problem for current time step
     * 
     * PV is solved independently (one-way coupling means p_l = 0 in PV eq.)
     */
    void solvePV();
    
    /**
     * @brief Solve PLC problem for current time step
     * 
     * Uses the polynomial BC coefficients from interpolatePVtoPLC()
     */
    void solvePLC();
    
    /**
     * @brief Transfer PV pressure to PLC boundary via polynomial fitting
     * 
     * Steps:
     * 1. Extract PV pressure along a vertical line (z-direction)
     * 2. Fit polynomial: p_v(z) = c0 + c1*z + c2*z^2 + ...
     * 3. Store coefficients in M_pvBCCoefficients
     * 4. Pass coefficients to PLC problem for BC evaluation
     */
    void interpolatePVtoPLC();
    
    /**
     * @brief Update solutions for next time step
     */
    void updateSolutions();
    
    // ========================================================================
    // Coefficient Access (for PLC BC evaluation)
    // ========================================================================
    
    /**
     * @brief Get polynomial coefficients for PV pressure BC
     * @return Const reference to coefficient vector [c0, c1, c2, ...]
     * 
     * The polynomial is: p_v(z) = c0 + c1*z + c2*z^2 + ... + cn*z^n
     * where z is normalized to [0,1] or uses actual z coordinates
     */
    const std::vector<scalar_type>& getPVBCCoefficients() const { 
        return M_pvBCCoefficients; 
    }
    
    /**
     * @brief Evaluate PV pressure at a given z coordinate
     * @param z Z-coordinate (vertical position along osteon axis)
     * @return Interpolated pressure value p_v(z)
     */
    scalar_type evaluatePVPressure(scalar_type z) const;
    
    /**
     * @brief Get the z-coordinate range for the interpolation
     * @param[out] z_min Minimum z
     * @param[out] z_max Maximum z
     */
    void getZRange(scalar_type& z_min, scalar_type& z_max) const;
    
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
     * @brief Export interpolation data for debugging/visualization
     * @param folder Output directory
     * @param frame Time frame number
     */
    void exportInterpolationData(const std::string& folder, int frame);
    
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
    inline TimeLoop* getTime() { return M_time; }
    
    /// Get leakage coefficient γ (for potential future two-way coupling)
    inline scalar_type getLeakageCoeff() const { return M_gamma; }
    
    /// Get coupling approach
    inline CouplingApproach getCouplingApproach() const { return M_couplingApproach; }
    
    /// Check if initialized
    inline bool isInitialized() const { return M_initialized; }

private:
    // ========================================================================
    // Private Helper Methods
    // ========================================================================
    
    /**
     * @brief Setup interpolation line profile based on geometry
     * 
     * Creates a vertical line through the domain center for z-sampling
     */
    void setupLineProfile();
    
    /**
     * @brief Update PLC boundary conditions with new coefficients
     * 
     * Passes polynomial coefficients to PLC problem's BC handler
     */
    void updatePLCBoundaryCoefficients();
    
    // ========================================================================
    // Member Variables - Domains
    // ========================================================================
    
    Bulk* M_bulkPV;           ///< PV (outer) bulk domain
    Bulk* M_bulkPLC;          ///< PLC (inner) bulk domain
    TimeLoop* M_time;         ///< Time manager
    
    // ========================================================================
    // Member Variables - Sub-problems
    // ========================================================================
    
    CoupledProblem* M_pvProblem;    ///< PV poroelasticity problem
    PLCProblem* M_plcProblem;       ///< PLC poroelasticity problem
    
    // ========================================================================
    // Member Variables - Interpolation
    // ========================================================================
    
    /// Interpolation manager (owned)
    std::unique_ptr<InterpolationManager> M_interpManager;
    
    /// Coupling approach
    CouplingApproach M_couplingApproach;
    
    /// Polynomial coefficients for PV pressure: p_v(z) = sum_i c_i * z^i
    std::vector<scalar_type> M_pvBCCoefficients;
    
    /// Z-coordinate range for normalization
    scalar_type M_z_min;
    scalar_type M_z_max;
    
    // ========================================================================
    // Member Variables - Coupling Parameters
    // ========================================================================
    
    /// Leakage coefficient γ (for potential two-way coupling)
    scalar_type M_gamma;
    
    // ========================================================================
    // Member Variables - Configuration
    // ========================================================================
    
    std::string M_section;           ///< Data file section
    bool M_initialized;              ///< Initialization flag
    size_type M_polynomialOrder;     ///< Order of polynomial fit
};

#endif // RUSSIANDOLLPROBLEM_H