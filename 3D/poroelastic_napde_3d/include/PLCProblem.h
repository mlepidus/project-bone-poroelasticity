// ============================================================================
// PLCProblem.h - Lacuno-canalicular porosity (PLC) poroelasticity problem
// ============================================================================
// Enhanced version with:
// - Polynomial BC from PV pressure (outer wall Dirichlet)
// - Coupling RHS term gamma * M * p_v for mass balance equation
// - Support for both one-way and two-way coupling
// ============================================================================
#ifndef PLCPROBLEM_H
#define PLCPROBLEM_H

#include "CoupledProblem.h"
#include <functional>

// Forward declaration
class RussianDollProblem;
class InterpolationManager;
struct PolynomialFit;
struct FieldPolynomials;

/**
 * @class PLCProblem
 * @brief Lacuno-canalicular porosity problem with inter-porosity coupling
 * 
 * Inherits from CoupledProblem and adds:
 * 1. Support for polynomial-based Dirichlet BC on outer wall (from PV pressure)
 * 2. Coupling RHS term: +gamma * M * p_v for the mass balance equation
 * 3. Access to polynomial coefficients for debugging/export
 * 
 * The governing equations are:
 * 
 *   K u_l^N - alpha_l C p_l^N = R_l^N                                      (momentum)
 *   alpha_l C^T u_dot_l^N + (H_l + gamma*M) p_l^N + (1/M_l) M p_dot_l^N 
 *       = Q_l^N + gamma*M*p_v^N                                            (mass)
 * 
 * Boundary Conditions:
 *   - Outer wall: p_l = p_v(z) where p_v(z) is polynomial from RussianDoll
 *   - Inner wall, top, bottom: as specified in data file
 * 
 * For the one-way coupling approach, the outer wall BC uses:
 *   p_l(x,y,z) = f(z) where f(z) = c0 + c1*z + c2*z^2 + ...
 * 
 * The coupling RHS adds +gamma * M * p_v where p_v is evaluated at each DOF
 * using the polynomial fit, allowing the leakage term to be included.
 */
class PLCProblem : public CoupledProblem {
public:
    // ========================================================================
    // Type Definitions
    // ========================================================================
    
    /**
     * @typedef BCCallback
     * @brief Function type for evaluating pressure BC as a function of z
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
     * @param section Data file section for PLC parameters (default: "bulkDataPLC/")
     */
    PLCProblem(const GetPot& dataFile, 
               Bulk* bulk = nullptr, 
               TimeLoop* time = nullptr,
               const std::string& section = "bulkDataPLC/");
    
    /// Virtual destructor
    virtual ~PLCProblem() = default;
    
    // ========================================================================
    // Boundary Condition Methods
    // ========================================================================
    
    /**
     * @brief Set callback for evaluating outer wall pressure BC
     * @param callback Function that evaluates p_v(z) at given z
     */
    void setOuterWallBCCallback(BCCallback callback);
    
    /**
     * @brief Check if outer wall BC callback is set
     */
    inline bool hasOuterWallBCCallback() const { return M_outerWallBCCallbackP != nullptr; }
    
    /**
     * @brief Set polynomial coefficients directly for outer wall BC
     * @param coefficients Polynomial coefficients [c0, c1, c2, ...]
     * @param z_min Minimum z-coordinate for normalization
     * @param z_max Maximum z-coordinate for normalization
     */
    void setOuterWallBCCoefficients(const std::vector<scalar_type>& coefficients,
                                    scalar_type z_min, scalar_type z_max);
    
    void setOuterWallDisplacementBCCoefficients(const std::vector<scalar_type>& coeffs_x,
        const std::vector<scalar_type>& coeffs_y,
        const std::vector<scalar_type>& coeffs_z,
        scalar_type z_min, scalar_type z_max);

    /**
     * @brief Get stored polynomial coefficients
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
    void computeDisplacementOnPLCDOFs();
    /**
     * @brief Get outer wall region ID
     */
    inline size_type getOuterWallRegion() const { return M_outerWallRegion; }
    
    /**
     * @brief Set outer wall region ID
     */
    void setOuterWallRegion(size_type region) { M_outerWallRegion = region; }
    
    /**
     * @brief Get z-range for BC evaluation
     */
    void getZRange(scalar_type& z_min, scalar_type& z_max) const {
        z_min = M_z_min;
        z_max = M_z_max;
    }
    
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
     * @brief Set the PV pressure vector (interpolated onto PLC mesh DOFs)
     * @param pv_on_PLC Pressure from PV problem evaluated at PLC DOF locations
     * 
     * This is used for the RHS coupling term: +gamma * M * p_v
     */
    void setCouplingSource(const scalarVectorPtr_Type& pv_on_PLC);
    
    /**
     * @brief Set PV pressure using polynomial coefficients
     * @param coefficients Polynomial coefficients [c0, c1, ...]
     * @param z_min Minimum z for normalization
     * @param z_max Maximum z for normalization
     * 
     * Evaluates polynomial at each PLC pressure DOF and stores the result
     */
    void setCouplingSourceFromPolynomial(const std::vector<scalar_type>& coefficients,
                                         scalar_type z_min, scalar_type z_max);
    
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
     * Builds standard poroelasticity matrices plus pressure mass matrix
     * for coupling term (+gamma * M * p_l - gamma * M * p_v)
     */
    virtual void assembleMatrix() override;
    
    /**
     * @brief Assemble global RHS
     * 
     * Calls base class assembleRHS(), adds coupling term (+gamma * M * p_v)
     */
    virtual void assembleRHS() override;
    
    /**
     * @brief Enforce boundary conditions with polynomial outer wall BC
     * @param firstTime If true, modify matrix and RHS; if false, only RHS
     */
    virtual void enforceStrongBC(bool firstTime) override;
    
    // ========================================================================
    // Additional Methods
    // ========================================================================
    
    /**
     * @brief Compute and add the coupling matrix term to system matrix
     * 
     * Adds +gamma * M to the (pressure, pressure) block
     */
    void addCouplingMatrix();
    
    /**
     * @brief Assemble Neumann BC contribution on outer wall
     * 
     * Adds -∫_Γ p_v * v·n dS to the velocity equation RHS
     */
    void assembleOuterWallNeumannRHS();
    

    /**
     * @brief Get pressure mass matrix (for coupling term)
     */
    inline sparseMatrixPtr_Type getPressureMassMatrix() const { return M_pressureMass; }
    
    /**
     * @brief Print coupling information for debugging
     */
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
     * @brief Enforce polynomial Dirichlet BC on outer wall
     * @param firstTime Modify matrix if true
     */
    void enforceOuterWallBC(bool firstTime);
    
        /**
     * @brief Enforce polynomial Dirichlet BC on outer wall
     * @param firstTime Modify matrix if true
     */
    void enforceOuterWallDisplacementBC(bool firstTime);

    /**
     * @brief Evaluate polynomial at a point
     * @param z Z-coordinate
     * @param coeffs Polynomial coefficients
     * @param z_min Minimum z for normalization
     * @param z_max Maximum z for normalization
     * @return Polynomial value
     */
    scalar_type evaluatePolynomial(scalar_type z,
                                   const std::vector<scalar_type>& coeffs,
                                   scalar_type z_min,
                                   scalar_type z_max) const;
    
    // ========================================================================
    // Member Variables
    // ========================================================================
    
    /// Leakage coefficient gamma [1/(Pa*s)]
    scalar_type M_gamma;
    
    /// Pressure from PV interpolated onto PLC mesh (for RHS coupling term)
    scalarVectorPtr_Type M_couplingSource;
    
    /// Pressure mass matrix for coupling term: M = integral N^T N dV
    sparseMatrixPtr_Type M_pressureMass;
    
    /// Flag: has pressure mass matrix been built?
    bool M_pressureMassBuilt;
    
    /// Flag: has coupling matrix been added to system?
    bool M_couplingMatrixAdded;
    
    /// Callback for evaluating outer wall BC: p_l = p_v(z)
    BCCallback M_outerWallBCCallbackP;

    /// Callback for evaluating outer wall displacement BC: 
    BCCallback M_outerWallBCCallbackX;
    BCCallback M_outerWallBCCallbackY;
    BCCallback M_outerWallBCCallbackZ;

    
    /// Polynomial coefficients for outer wall BC
    std::vector<scalar_type> M_bcCoefficients; //pressure coefficients
    
    std::vector<scalar_type> M_dispXCoefficients;  // Displacement X coefficients
    std::vector<scalar_type> M_dispYCoefficients;  // Displacement Y coefficients
    std::vector<scalar_type> M_dispZCoefficients;  // Displacement Z coefficients
    /// Z-range for polynomial normalization
    scalar_type M_z_min;
    scalar_type M_z_max;
    
    /// Outer wall region ID (typically region 0 = outer surface)
    size_type M_outerWallRegion;
    
    /// Flag: use polynomial BC on outer wall
    bool M_usePolynomialBC;
    bool M_useDisplacementBC;

};

#endif // PLCPROBLEM_H