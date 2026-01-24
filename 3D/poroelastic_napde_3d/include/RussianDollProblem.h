// ============================================================================
// RussianDollProblem.h - Dual-porosity (PV + PLC) coupled problem
// ============================================================================
// Enhanced version with:
// - Line interpolation for pressure and displacement
// - Polynomial coefficient computation and storage
// - One-way and two-way coupling support
// - Complete workflow: solve PV -> interpolate -> modify RHS -> solve PLC
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
 * This class orchestrates the coupling workflow:
 * 
 * For ONE_WAY coupling:
 * 1. Solve PV problem (uncoupled from PLC, i.e., p_l = 0 in PV equation)
 * 2. Extract PV pressure and displacement along a vertical line (z-direction)
 * 3. Fit polynomials: p_v(z) = c0 + c1*z + c2*z^2 + ..., similarly for u_v(z)
 * 4. Store polynomial coefficients for BC evaluation
 * 5. Update PLC outer wall BC: p_l|_outer = p_v(z)
 * 6. Add coupling RHS: +gamma * M * p_v
 * 7. Solve PLC problem
 * 8. Export results
 * 
 * The governing equations are:
 * 
 * PV (Vascular Porosity):
 *   K u_v^N - alpha_v C p_v^N = R_v^N
 *   alpha_v C^T u_dot_v^N + H_v p_v^N + (1/M_v) M p_dot_v^N = Q_v^N
 * 
 * PLC (Lacuno-canalicular Porosity):
 *   K u_l^N - alpha_l C p_l^N = R_l^N
 *   alpha_l C^T u_dot_l^N + (H_l + gamma*M) p_l^N + (1/M_l) M p_dot_l^N 
 *       = Q_l^N + gamma*M*p_v^N
 *   with BC: p_l = p_v(z) on outer wall (Dirichlet from interpolated PV pressure)
 * 
 * Key Features:
 * - Polynomial coefficients from PV -> PLC transfer are explicitly stored
 * - Supports pressure and displacement interpolation
 * - PLCProblem's BC class uses coefficients through callbacks
 * - Complete export of interpolation data for debugging
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
     */
    void setPLCProblem(PLCProblem* plcProblem);
    
    /**
     * @brief Initialize both problems and interpolation structures
     */
    void initialize();
    
    // ========================================================================
    // Solution Methods - Staggered One-Way Coupling
    // ========================================================================
    
    /**
     * @brief Perform one time step with one-way coupling (PV -> PLC)
     * 
     * Algorithm:
     * 1. Solve PV problem (standalone, p_l = 0)
     * 2. Extract PV pressure along vertical line
     * 3. Fit polynomial p_v(z)
     * 4. Update PLC boundary coefficients
     * 5. Add coupling RHS: +gamma * M * p_v
     * 6. Solve PLC problem with polynomial BC on outer wall
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
     * @brief Transfer PV solution to PLC boundary via polynomial fitting
     * 
     * Steps:
     * 1. Extract PV pressure along a vertical line (z-direction)
     * 2. Extract PV displacement along the same line
     * 3. Fit polynomials: p_v(z), u_v_x(z), u_v_y(z), u_v_z(z)
     * 4. Store coefficients for BC evaluation
     * 5. Pass coefficients to PLC problem
     */
    void interpolatePVtoPLC();
    
    /**
     * @brief Update solutions for next time step
     */
    void updateSolutions();
    
    // ========================================================================
    // Coefficient Access
    // ========================================================================
    
    /**
     * @brief Get polynomial coefficients for PV pressure BC
     * @return Const reference to coefficient vector [c0, c1, c2, ...]
     */
    const std::vector<scalar_type>& getPressureCoefficients() const { 
        return M_pressureCoefficients; 
    }
    
    /**
     * @brief Get polynomial coefficients for PV displacement
     * @param component 0=x, 1=y, 2=z
     * @return Const reference to coefficient vector
     */
    const std::vector<scalar_type>& getDisplacementCoefficients(int component) const;
    
    /**
     * @brief Evaluate PV pressure at a given z coordinate
     * @param z Z-coordinate (vertical position along osteon axis)
     * @return Interpolated pressure value p_v(z)
     */
    scalar_type evaluatePVPressure(scalar_type z) const;
    
    /**
     * @brief Evaluate PV displacement at a given z coordinate
     * @param z Z-coordinate
     * @param[out] ux X-component
     * @param[out] uy Y-component
     * @param[out] uz Z-component
     */
    void evaluatePVDisplacement(scalar_type z, scalar_type& ux, scalar_type& uy, scalar_type& uz) const;
    
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
     * 
     * Outputs:
     * - Polynomial coefficients
     * - Sampled values along the interpolation line
     * - BC values at sample z-coordinates
     */
    void exportInterpolationData(const std::string& folder, int frame);
    
    /**
     * @brief Export polynomial coefficients to file
     * @param filename Output file path
     */
    void exportCoefficients(const std::string& filename) const;
    
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
    
    /// Get leakage coefficient gamma
    inline scalar_type getLeakageCoeff() const { return M_gamma; }
    
    /// Get coupling approach
    inline CouplingApproach getCouplingApproach() const { return M_couplingApproach; }
    
    /// Check if initialized
    inline bool isInitialized() const { return M_initialized; }
    
    /// Get polynomial order
    inline size_type getPolynomialOrder() const { return M_polynomialOrder; }

private:
    // ========================================================================
    // Private Helper Methods
    // ========================================================================
    
    /**
     * @brief Setup interpolation line profile based on geometry
     */
    void setupLineProfile();
    
    /**
     * @brief Update PLC boundary conditions with new coefficients
     */
    void updatePLCBoundaryCoefficients();
    
    /**
     * @brief Compute PV pressure at all PLC pressure DOF locations
     * 
     * Uses polynomial coefficients to evaluate p_v(z) at each PLC DOF
     */
    void computePVPressureOnPLCDOFs();
    
    /**
     * @brief Print polynomial coefficient information
     */
    void printCoefficients() const;
    
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
    std::vector<scalar_type> M_pressureCoefficients;
    
    /// Polynomial coefficients for PV displacement components
    std::vector<scalar_type> M_displacementXCoefficients;
    std::vector<scalar_type> M_displacementYCoefficients;
    std::vector<scalar_type> M_displacementZCoefficients;
    
    /// Z-coordinate range for normalization
    scalar_type M_z_min;
    scalar_type M_z_max;
    
    // ========================================================================
    // Member Variables - Coupling Parameters
    // ========================================================================
    
    /// Leakage coefficient gamma
    scalar_type M_gamma;
    
    // ========================================================================
    // Member Variables - Configuration
    // ========================================================================
    
    std::string M_section;           ///< Data file section
    bool M_initialized;              ///< Initialization flag
    size_type M_polynomialOrder;     ///< Order of polynomial fit
    
    /// Store PV pressure evaluated at PLC DOFs for coupling RHS
    scalarVectorPtr_Type M_pvPressureOnPLC;
};

#endif // RUSSIANDOLLPROBLEM_H