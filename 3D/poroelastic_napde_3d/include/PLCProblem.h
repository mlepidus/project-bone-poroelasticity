// ============================================================================
// PLCProblem.h - Lacuno-canalicular porosity (PLC) poroelasticity problem
// ============================================================================
// This version supports polynomial-based boundary conditions from RussianDoll
// coupling. The outer wall BC can be set via polynomial coefficients p_v(z).
// ============================================================================
#ifndef PLCPROBLEM_H
#define PLCPROBLEM_H

#include "CoupledProblem.h"
#include <functional>

// Forward declaration
class RussianDollProblem;

/**
 * @class PLCProblem
 * @brief Lacuno-canalicular porosity problem with inter-porosity coupling
 * 
 * Inherits from CoupledProblem and adds:
 * 1. Support for polynomial-based Dirichlet BC on outer wall (from PV pressure)
 * 2. Optional coupling source term (for two-way coupling, not used in one-way)
 * 
 * The governing equations are:
 * 
 *   K u_l^N - α_l C p_l^N = R_l^N                                     (momentum)
 *   α_l C^T u̇_l^N + H_l p_l^N + (1/M_l) M ṗ_l^N = Q_l^N              (mass)
 * 
 * Boundary Conditions:
 *   - Outer wall: p_l = p_v(z) where p_v(z) is polynomial from RussianDoll
 *   - Inner wall, top, bottom: as specified in data file
 * 
 * For the one-way coupling approach, the outer wall BC uses:
 *   p_l(x,y,z) = f(z) where f(z) = c0 + c1*z + c2*z^2 + ...
 * 
 * Usage:
 * @code
 *   PLCProblem plc(dataFile, bulkPLC, time, "plc/");
 *   plc.addElastPB(&elastPLC);
 *   plc.addDarcyPB(&darcyPLC);
 *   
 *   // Set polynomial BC callback from RussianDoll
 *   plc.setOuterWallBCCallback([&russian](scalar_type z) {
 *       return russian.evaluatePVPressure(z);
 *   });
 *   
 *   // In time loop:
 *   plc.assembleRHS();
 *   plc.enforceStrongBC(false);  // Uses polynomial callback for outer wall
 *   plc.solve();
 * @endcode
 */
class PLCProblem : public CoupledProblem {
public:
    // ========================================================================
    // Type Definitions
    // ========================================================================
    
    /**
     * @typedef BCCallback
     * @brief Function type for evaluating pressure BC as a function of z
     * 
     * The callback takes z-coordinate and returns pressure value.
     * This allows RussianDoll to provide polynomial evaluation.
     */
    using BCCallback = std::function<scalar_type(scalar_type z)>;
    
    // ========================================================================
    // Construction / Destruction
    // ========================================================================
    
    /**
     * @brief Construct PLC problem from data file
     * @param dataFile GetPot parameter file
     * @param bulk Pointer to PLC bulk domain (inner mesh)
     * @param time Pointer to time loop manager
     * @param section Data file section for PLC parameters (default: "plc/")
     */
    PLCProblem(const GetPot& dataFile, 
               Bulk* bulk = nullptr, 
               TimeLoop* time = nullptr,
               const std::string& section = "plc/");
    
    /// Virtual destructor
    virtual ~PLCProblem() = default;
    
    // ========================================================================
    // Boundary Condition Methods
    // ========================================================================
    
    /**
     * @brief Set callback for evaluating outer wall pressure BC
     * @param callback Function that evaluates p_v(z) at given z
     * 
     * This callback is used in enforceStrongBC() to set Dirichlet values
     * on the outer wall based on the polynomial fit from PV pressure.
     * 
     * Example:
     * @code
     *   plc.setOuterWallBCCallback([coeffs, z_min, z_max](scalar_type z) {
     *       scalar_type t = (z - z_min) / (z_max - z_min);  // normalize
     *       scalar_type val = coeffs[0];
     *       scalar_type t_pow = t;
     *       for (size_t i = 1; i < coeffs.size(); ++i) {
     *           val += coeffs[i] * t_pow;
     *           t_pow *= t;
     *       }
     *       return val;
     *   });
     * @endcode
     */
    void setOuterWallBCCallback(BCCallback callback);
    
    /**
     * @brief Check if outer wall BC callback is set
     * @return true if callback is set and valid
     */
    inline bool hasOuterWallBCCallback() const { return M_outerWallBCCallback != nullptr; }
    
    /**
     * @brief Set polynomial coefficients directly for outer wall BC
     * @param coefficients Polynomial coefficients [c0, c1, c2, ...]
     * @param z_min Minimum z-coordinate for normalization
     * @param z_max Maximum z-coordinate for normalization
     * 
     * Alternative to setOuterWallBCCallback() - creates internal callback
     */
    void setOuterWallBCCoefficients(const std::vector<scalar_type>& coefficients,
                                    scalar_type z_min, scalar_type z_max);
    
    /**
     * @brief Get stored polynomial coefficients
     * @return Const reference to coefficient vector
     */
    const std::vector<scalar_type>& getOuterWallBCCoefficients() const {
        return M_bcCoefficients;
    }
    
    /**
     * @brief Evaluate outer wall BC at a point
     * @param z Z-coordinate
     * @return Pressure value from callback or coefficients
     */
    scalar_type evaluateOuterWallBC(scalar_type z) const;
    
    /**
     * @brief Get outer wall region ID
     * @return Region ID for outer wall boundary
     */
    inline size_type getOuterWallRegion() const { return M_outerWallRegion; }
    
    /**
     * @brief Set outer wall region ID
     * @param region Region ID
     */
    void setOuterWallRegion(size_type region) { M_outerWallRegion = region; }
    
    // ========================================================================
    // Coupling Methods (for optional two-way coupling)
    // ========================================================================
    
    /**
     * @brief Set leakage coefficient γ for inter-porosity coupling
     * @param gamma Leakage coefficient [1/(Pa·s)]
     * 
     * This is for potential two-way coupling (not used in one-way mode)
     */
    void setLeakageCoefficient(scalar_type gamma);
    
    /// Get leakage coefficient
    inline scalar_type getLeakageCoefficient() const { return M_gamma; }
    
    /**
     * @brief Set coupling source term (p_v interpolated onto PLC mesh)
     * @param pv_on_PLC Pressure from PV problem interpolated to PLC DOFs
     * 
     * For two-way coupling: adds +γ * M * p_v to RHS
     * Not used in one-way coupling mode
     */
    void setCouplingSource(const scalarVectorPtr_Type& pv_on_PLC);
    
    /// Get coupling source
    inline scalarVectorPtr_Type getCouplingSource() const { return M_couplingSource; }
    
    /// Check if coupling source is set
    inline bool hasCouplingSource() const { 
        return M_couplingSource != nullptr && M_couplingSource->size() > 0; 
    }
    
    // ========================================================================
    // Overridden Assembly Methods
    // ========================================================================
    
    /**
     * @brief Assemble global system matrix
     * 
     * Same as base class, optionally builds pressure mass matrix for coupling
     */
    virtual void assembleMatrix() override;
    
    /**
     * @brief Assemble global RHS
     * 
     * Calls base class assembleRHS(), optionally adds coupling term
     */
    virtual void assembleRHS() override;
    
    /**
     * @brief Enforce boundary conditions with polynomial outer wall BC
     * @param firstTime If true, modify matrix and RHS; if false, only RHS
     * 
     * Uses the outer wall BC callback to evaluate pressure at each DOF
     */
    virtual void enforceStrongBC(bool firstTime) override;
    
    // ========================================================================
    // Additional Getters
    // ========================================================================
    
    /// Get pressure mass matrix (for coupling term)
    inline sparseMatrixPtr_Type getPressureMassMatrix() const { return M_pressureMass; }

protected:
    // ========================================================================
    // Protected Methods
    // ========================================================================
    
    /**
     * @brief Build pressure mass matrix for coupling term
     * 
     * Computes: M = ∫ N^T N dV where N are pressure shape functions
     */
    void buildPressureMassMatrix();
    
    /**
     * @brief Assemble the inter-porosity coupling RHS term
     * 
     * Computes: +γ * M * p_v and adds to RHS (for two-way coupling)
     */
    void assembleCouplingRHS();
    
    /**
     * @brief Enforce polynomial Dirichlet BC on outer wall
     * @param firstTime Modify matrix if true
     * 
     * Uses M_outerWallBCCallback to evaluate p_v(z) at each DOF on outer wall
     */
    void enforceOuterWallBC(bool firstTime);
    
    // ========================================================================
    // Member Variables
    // ========================================================================
    
    /// Leakage coefficient γ [1/(Pa·s)] (for two-way coupling)
    scalar_type M_gamma;
    
    /// Pressure from PV interpolated onto PLC mesh (for two-way coupling)
    scalarVectorPtr_Type M_couplingSource;
    
    /// Pressure mass matrix for coupling term: M = ∫ N^T N dV
    sparseMatrixPtr_Type M_pressureMass;
    
    /// Flag: has pressure mass matrix been built?
    bool M_pressureMassBuilt;
    
    /// Callback for evaluating outer wall BC: p_l = p_v(z)
    BCCallback M_outerWallBCCallback;
    
    /// Polynomial coefficients for outer wall BC (alternative to callback)
    std::vector<scalar_type> M_bcCoefficients;
    
    /// Z-range for polynomial normalization
    scalar_type M_z_min;
    scalar_type M_z_max;
    
    /// Outer wall region ID (typically region 0 = outer surface)
    size_type M_outerWallRegion;
    
    /// Flag: use polynomial BC on outer wall
    bool M_usePolynomialBC;
};

#endif // PLCPROBLEM_H