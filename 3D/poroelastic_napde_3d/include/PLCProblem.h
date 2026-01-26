// ============================================================================
// PLCProblem.h - Lacuno-canalicular porosity (PLC) poroelasticity problem
// ============================================================================
// Simplified version that delegates polynomial BC handling to the BC class.
// 
// Key changes from previous version:
// - Polynomial pressure/displacement BCs are set via BC class callbacks
// - Uses standard enforceStrongBC() from ElastProblem (Nitsche or strong)
// - Uses standard naturalRHS() from DarcyOperatorsBD for pressure Neumann
// - Only adds coupling term (+gamma * M * p_v) to RHS
// ============================================================================
#ifndef PLCPROBLEM_H
#define PLCPROBLEM_H

#include "CoupledProblem.h"
#include <functional>

// Forward declarations
class RussianDollProblem;
class InterpolationManager;

/**
 * @class PLCProblem
 * @brief Lacuno-canalicular porosity problem with inter-porosity coupling
 * 
 * Inherits from CoupledProblem and adds:
 * 1. Support for polynomial-based BCs via BC class callbacks
 * 2. Coupling RHS term: +gamma * M * p_v for the mass balance equation
 * 3. Coupling matrix term: +gamma * M added to pressure block
 * 
 * The governing equations are:
 * 
 *   K u_l^N - alpha_l C p_l^N = R_l^N                                      (momentum)
 *   alpha_l C^T u_dot_l^N + (H_l + gamma*M) p_l^N + (1/M_l) M p_dot_l^N 
 *       = Q_l^N + gamma*M*p_v^N                                            (mass)
 * 
 * Boundary Conditions (handled by BC class with callbacks):
 *   - Outer wall displacement: Dirichlet u_l = u_v(z) via BCDiriVec callback
 *   - Outer wall pressure: Neumann p_l = p_v(z) via BCNeum callback  
 *   - Other boundaries: as specified in data file bcflag
 */
class PLCProblem : public CoupledProblem {
public:
    // ========================================================================
    // Type Definitions
    // ========================================================================
    
    using ScalarCallback = std::function<scalar_type(scalar_type z)>;
    using VectorCallback = std::function<bgeot::base_node(scalar_type z)>;
    
    // ========================================================================
    // Construction / Destruction
    // ========================================================================
    
    /**
     * @brief Construct PLC problem from data file
     * @param dataFile GetPot parameter file
     * @param bulk Pointer to PLC bulk domain (inner mesh)
     * @param time Pointer to time loop manager
     * @param section Data file section for PLC parameters (default: "bulkDataPLC/")
     */
    PLCProblem(const GetPot& dataFile, 
               Bulk* bulk = nullptr, 
               TimeLoop* time = nullptr,
               const std::string& section = "bulkDataPLC/");
    
    /// Virtual destructor
    virtual ~PLCProblem() = default;
    
    // ========================================================================
    // Polynomial BC Setup (delegates to BC class)
    // ========================================================================
    
    /**
     * @brief Set polynomial pressure BC on outer wall
     * @param coefficients Polynomial coefficients [c0, c1, c2, ...]
     * @param z_min Minimum z-coordinate for normalization
     * @param z_max Maximum z-coordinate for normalization
     * 
     * This sets up the BC class callback for pressure (Neumann).
     * The bcflag for outer wall (region 0) should be 1 (Neumann) for Darcy.
     */
    void setOuterWallPressureBC(const std::vector<scalar_type>& coefficients,
                                scalar_type z_min, scalar_type z_max);
    
    /**
     * @brief Set polynomial displacement BC on outer wall
     * @param coeffs_x X-component polynomial coefficients
     * @param coeffs_y Y-component polynomial coefficients  
     * @param coeffs_z Z-component polynomial coefficients
     * @param z_min Minimum z-coordinate for normalization
     * @param z_max Maximum z-coordinate for normalization
     * 
     * This sets up the BC class callback for displacement (Dirichlet).
     * The bcflag for outer wall (region 0) should be 0 (Dirichlet) for mecc.
     */
    void setOuterWallDisplacementBC(const std::vector<scalar_type>& coeffs_x,
                                    const std::vector<scalar_type>& coeffs_y,
                                    const std::vector<scalar_type>& coeffs_z,
                                    scalar_type z_min, scalar_type z_max);
    
    /**
     * @brief Clear polynomial BCs on outer wall (revert to data file values)
     */
    void clearOuterWallPolynomialBCs();
    
    // ========================================================================
    // Coupling Methods
    // ========================================================================
    
    /**
     * @brief Set leakage coefficient gamma for inter-porosity coupling
     * @param gamma Leakage coefficient [1/(Pa*s)]
     */
    void setLeakageCoefficient(scalar_type gamma);
    
    /// Get leakage coefficient
    inline scalar_type getLeakageCoefficient() const { return M_gamma; }
    
    /**
     * @brief Set the PV pressure vector for coupling RHS
     * @param pv_on_PLC Pressure from PV evaluated at PLC DOF locations
     * 
     * This is used for the RHS coupling term: +gamma * M * p_v
     */
    void setCouplingPressure(const scalarVectorPtr_Type& pv_on_PLC);
    
    /**
     * @brief Set PV pressure using polynomial coefficients
     * @param coefficients Polynomial coefficients [c0, c1, ...]
     * @param z_min Minimum z for normalization
     * @param z_max Maximum z for normalization
     * 
     * Evaluates polynomial at each PLC pressure DOF and stores the result.
     */
    void setCouplingPressureFromPolynomial(const std::vector<scalar_type>& coefficients,
                                           scalar_type z_min, scalar_type z_max);
    
    /// Check if coupling source is set
    inline bool hasCouplingPressure() const { 
        return M_couplingPressure != nullptr && M_couplingPressure->size() > 0; 
    }
    
    // ========================================================================
    // Overridden Assembly Methods
    // ========================================================================
    
    /**
     * @brief Assemble global system matrix
     * 
     * Calls base class, then builds pressure mass matrix for coupling.
     */
    virtual void assembleMatrix() override;
    
    /**
     * @brief Assemble global RHS
     * 
     * Calls base class assembleRHS(), then adds coupling term (+gamma * M * p_v).
     */
    virtual void assembleRHS() override;
    
    /**
     * @brief Enforce boundary conditions
     * @param firstTime If true, modify matrix and RHS; if false, only RHS
     * 
     * Calls base class which uses standard Nitsche/strong BC enforcement.
     * The polynomial callbacks are already set in the BC class.
     * Also adds coupling matrix term (+gamma * M) on first call.
     */
    virtual void enforceStrongBC(bool firstTime) override;
    
    // ========================================================================
    // Getters
    // ========================================================================
    
    /// Get outer wall region ID
    inline size_type getOuterWallRegion() const { return M_outerWallRegion; }
    
    /// Set outer wall region ID
    void setOuterWallRegion(size_type region) { M_outerWallRegion = region; }
    
    /// Get z-range for BC evaluation
    void getZRange(scalar_type& z_min, scalar_type& z_max) const {
        z_min = M_z_min;
        z_max = M_z_max;
    }
    
    /// Get pressure mass matrix (for coupling term)
    inline sparseMatrixPtr_Type getPressureMassMatrix() const { return M_pressureMass; }
    
    /// Print coupling information for debugging
    void printCouplingInfo() const;

protected:
    // ========================================================================
    // Protected Methods
    // ========================================================================
    
    /**
     * @brief Build pressure mass matrix for coupling term
     * 
     * Computes: M = integral N^T N dV where N are pressure shape functions
     */
    void buildPressureMassMatrix();
    
    /**
     * @brief Assemble the inter-porosity coupling RHS term
     * 
     * Computes: +gamma * M * p_v and adds to RHS
     */
    void assembleCouplingRHS();
    
    /**
     * @brief Add coupling matrix term to system
     * 
     * Adds +gamma * M to the (pressure, pressure) block
     */
    void addCouplingMatrix();
    
    // ========================================================================
    // Member Variables
    // ========================================================================
    
    /// Leakage coefficient gamma [1/(Pa*s)]
    scalar_type M_gamma;
    
    /// Pressure from PV interpolated onto PLC mesh (for RHS coupling term)
    scalarVectorPtr_Type M_couplingPressure;
    
    /// Pressure mass matrix for coupling term: M = integral N^T N dV
    sparseMatrixPtr_Type M_pressureMass;
    
    /// Flag: has pressure mass matrix been built?
    bool M_pressureMassBuilt;
    
    /// Flag: has coupling matrix been added to system?
    bool M_couplingMatrixAdded;
    
    /// Z-range for polynomial BC normalization
    scalar_type M_z_min;
    scalar_type M_z_max;
    
    /// Outer wall region ID (typically region 0 = outer surface)
    size_type M_outerWallRegion;
    
    /// Store polynomial coefficients for debugging/export
    std::vector<scalar_type> M_pressureCoefficients;
    std::vector<scalar_type> M_dispXCoefficients;
    std::vector<scalar_type> M_dispYCoefficients;
    std::vector<scalar_type> M_dispZCoefficients;
};

#endif // PLCPROBLEM_H