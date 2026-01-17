// ============================================================================
// PLCProblem.h - Lacuno-canalicular porosity (PLC) poroelasticity problem
// ============================================================================
#ifndef PLCPROBLEM_H
#define PLCPROBLEM_H

#include "CoupledProblem.h"

/**
 * @class PLCProblem
 * @brief Lacuno-canalicular porosity problem with inter-porosity coupling
 * 
 * Inherits from CoupledProblem and adds the coupling term from vascular
 * porosity (PV). The governing equations are:
 * 
 *   K u_l^N - α_l C p_l^N = R_l^N                                     (momentum)
 *   α_l C^T u̇_l^N + (H_l + γM) p_l^N + (1/M_l) M ṗ_l^N = Q_l^N + γM p_v^N  (mass)
 * 
 * The additional term +γM p_v^N represents fluid exchange from the vascular
 * network and is set via setCouplingSource().
 * 
 * Usage:
 * @code
 *   PLCProblem plc(dataFile, bulkPLC, time, "plc/");
 *   plc.addElastPB(&elastPLC);
 *   plc.addDarcyPB(&darcyPLC);
 *   plc.setLeakageCoefficient(gamma);
 *   
 *   // In time loop:
 *   plc.setCouplingSource(pv_interpolated_to_plc);
 *   plc.assembleRHS();  // Includes coupling term
 *   plc.solve();
 * @endcode
 */
class PLCProblem : public CoupledProblem {
public:
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
    // Coupling Methods (Inter-porosity exchange with PV)
    // ========================================================================
    
    /**
     * @brief Set leakage coefficient γ for inter-porosity coupling
     * @param gamma Leakage coefficient [1/(Pa·s)]
     * 
     * This controls the strength of fluid exchange between PV and PLC networks.
     * The coupling term is: γ * (p_v - p_l)
     */
    void setLeakageCoefficient(scalar_type gamma);
    
    /**
     * @brief Get leakage coefficient
     * @return Current leakage coefficient
     */
    inline scalar_type getLeakageCoefficient() const { return M_gamma; }
    
    /**
     * @brief Set coupling source term (p_v interpolated onto PLC mesh)
     * @param pv_on_PLC Pressure from PV problem interpolated to PLC DOFs
     * 
     * This vector is used in assembleRHS() to add the term: +γ * M * p_v
     * Must be called before assembleRHS() at each time step.
     */
    void setCouplingSource(const scalarVectorPtr_Type& pv_on_PLC);
    
    /**
     * @brief Get the current coupling source
     * @return Pointer to coupling source vector (p_v on PLC mesh)
     */
    inline scalarVectorPtr_Type getCouplingSource() const { return M_couplingSource; }
    
    /**
     * @brief Check if coupling source has been set
     * @return true if coupling source is available
     */
    inline bool hasCouplingSource() const { 
        return M_couplingSource != nullptr && M_couplingSource->size() > 0; 
    }
    
    // ========================================================================
    // Overridden Assembly Methods
    // ========================================================================
    
    /**
     * @brief Assemble global system matrix
     * 
     * Same as base class, but also builds pressure mass matrix for coupling
     */
    virtual void assembleMatrix() override;
    
    /**
     * @brief Assemble global RHS including coupling term
     * 
     * Calls base class assembleRHS(), then adds: +γ * M * p_v
     * to the pressure equation (mass balance) RHS.
     */
    virtual void assembleRHS() override;
    
    // ========================================================================
    // Additional Getters
    // ========================================================================
    
    /// Get pressure mass matrix (used for coupling term)
    inline sparseMatrixPtr_Type getPressureMassMatrix() const { return M_pressureMass; }

protected:
    // ========================================================================
    // Protected Methods
    // ========================================================================
    
    /**
     * @brief Build pressure mass matrix for coupling term
     * 
     * Computes: M = ∫ N^T N dV where N are pressure shape functions
     * This is used for the leakage term: γ * M * p_v
     */
    void buildPressureMassMatrix();
    
    /**
     * @brief Assemble the inter-porosity coupling RHS term
     * 
     * Computes: +γ * M * p_v and adds to the Darcy sub-system RHS
     * at the pressure DOF locations.
     */
    void assembleCouplingRHS();
    
    // ========================================================================
    // Member Variables
    // ========================================================================
    
    /// Leakage coefficient γ [1/(Pa·s)]
    scalar_type M_gamma;
    
    /// Pressure from PV interpolated onto PLC mesh
    scalarVectorPtr_Type M_couplingSource;
    
    /// Pressure mass matrix for coupling term: M = ∫ N^T N dV
    sparseMatrixPtr_Type M_pressureMass;
    
    /// Flag: has pressure mass matrix been built?
    bool M_pressureMassBuilt;
};

#endif // PLCPROBLEM_H