// ============================================================================
// RussianDollProblem.h - Dual-porosity (PV + PLC) coupled problem
// ============================================================================
// Simplified and corrected version for one-way coupling:
// 
// Workflow at each time step:
// 1. Solve PV problem (standalone, p_l = 0 in coupling term)
// 2. Extract PV pressure and displacement along vertical line
// 3. Fit polynomials: p_v(z), u_v(z)
// 4. Update PLC boundary conditions via BC class callbacks
// 5. Set coupling RHS term: +gamma * M * p_v
// 6. Solve PLC problem
// 7. Export results
// ============================================================================
#ifndef RUSSIANDOLLPROBLEM_H
#define RUSSIANDOLLPROBLEM_H

#include "Core.h"
#include "LinearSystem.h"
#include "CoupledProblem.h"
#include "PLCProblem.h"
#include "InterpolationManager.h"
#include "TimeLoop.h"
#include <memory>

/**
 * @class RussianDollProblem
 * @brief Manages dual-porosity coupling between PV (vascular) and PLC (lacuno-canalicular)
 * 
 * One-way coupling: PV -> PLC
 * - PV is solved independently (p_l = 0)
 * - PV solution is interpolated along a line and fit to polynomials
 * - Polynomials are used as BCs for PLC outer wall
 * - Coupling term +gamma*M*p_v is added to PLC RHS
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
     * @brief Initialize interpolation structures
     * Must be called after problems are registered
     */
    void initialize();
    
    // ========================================================================
    // Core Coupling Methods
    // ========================================================================
    
    /**
     * @brief Interpolate PV solution and update PLC boundary conditions
     * 
     * This is the key coupling method:
     * 1. Get PV pressure and displacement solutions
     * 2. Extract values along the interpolation line (z-axis)
     * 3. Fit polynomials to extracted data
     * 4. Pass polynomial coefficients to PLC problem
     * 5. PLC uses these via BC class callbacks
     */
    void interpolatePVtoPLC();
    
    /**
     * @brief Update solutions for next time step
     */
    void updateSolutions();
    
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
     * @brief Export interpolation data for debugging
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
    
    /**
     * @brief Print coupling information for debugging
     */
    void printCouplingInfo() const;
    
    // ========================================================================
    // Getters
    // ========================================================================
    
    inline CoupledProblem* getPVProblem() { return M_pvProblem; }
    inline PLCProblem* getPLCProblem() { return M_plcProblem; }
    inline InterpolationManager* getInterpolationManager() { return M_interpManager.get(); }
    inline TimeLoop* getTime() { return M_time; }
    inline Bulk* getBulkPV() { return M_bulkPV; }
    inline Bulk* getBulkPLC() { return M_bulkPLC; }
    
    /// Get polynomial coefficients for pressure
    const std::vector<scalar_type>& getPressureCoefficients() const { 
        return M_pressureCoefficients; 
    }
    
    /// Get polynomial coefficients for displacement component
    const std::vector<scalar_type>& getDisplacementCoefficients(int component) const;
    
    /// Get z-range
    void getZRange(scalar_type& z_min, scalar_type& z_max) const {
        z_min = M_z_min;
        z_max = M_z_max;
    }
    
    /// Check if initialized
    inline bool isInitialized() const { return M_initialized; }

private:
    // ========================================================================
    // Private Helper Methods
    // ========================================================================
    
    /**
     * @brief Evaluate polynomial at normalized coordinate
     * @param coeffs Polynomial coefficients
     * @param z Physical z-coordinate
     * @return Polynomial value
     */
    scalar_type evaluatePolynomial(const std::vector<scalar_type>& coeffs, 
                                   scalar_type z) const;
    
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
    
    /// Polynomial coefficients for PV pressure
    std::vector<scalar_type> M_pressureCoefficients;
    
    /// Polynomial coefficients for PV displacement components
    std::vector<scalar_type> M_displacementXCoefficients;
    std::vector<scalar_type> M_displacementYCoefficients;
    std::vector<scalar_type> M_displacementZCoefficients;
    
    /// Z-coordinate range for normalization
    scalar_type M_z_min;
    scalar_type M_z_max;
    
    /// Polynomial order
    size_type M_polynomialOrder;
    
    // ========================================================================
    // Member Variables - Configuration
    // ========================================================================
    
    std::string M_section;           ///< Data file section
    bool M_initialized;              ///< Initialization flag
    
    /// Outer wall region ID for PLC
    size_type M_outerWallRegion;
};

#endif // RUSSIANDOLLPROBLEM_H