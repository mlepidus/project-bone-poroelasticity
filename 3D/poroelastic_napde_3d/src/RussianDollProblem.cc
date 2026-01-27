// ============================================================================
// RussianDollProblem.cc - Implementation of dual-porosity coupling
// ============================================================================

#include "../include/RussianDollProblem.h"
#include <iostream>
#include <fstream>
#include <iomanip>

// ============================================================================
// Constructor
// ============================================================================
RussianDollProblem::RussianDollProblem(const GetPot& dataFile, 
                                       Bulk* bulkPV, 
                                       Bulk* bulkPLC,
                                       TimeLoop* time)
    : M_bulkPV(bulkPV),
      M_bulkPLC(bulkPLC),
      M_time(time),
      M_pvProblem(nullptr),
      M_plcProblem(nullptr),
      M_interpManager(nullptr),
      M_z_min(0.0),
      M_z_max(1.0),
      M_polynomialOrder(5),
      M_section("russian_doll/"),
      M_initialized(false),
      M_outerWallRegion(0)
{
    // Read configuration from data file
    M_polynomialOrder = dataFile((M_section + "interpolation/line0/polynomial_order").c_str(), 5);
    M_outerWallRegion = dataFile("bulkDataPLC/outer_wall_region", 0);
    
    // Create interpolation manager
    M_interpManager = std::make_unique<InterpolationManager>(dataFile, bulkPV, bulkPLC);
    #ifdef VERBOSE
    std::cout << "=====================================================" << std::endl;
    std::cout << " RussianDollProblem Created" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Polynomial order: " << M_polynomialOrder << std::endl;
    std::cout << "  Outer wall region for PLC: " << M_outerWallRegion << std::endl;
    std::cout << "=====================================================" << std::endl;
    #endif

}

// ============================================================================
// Problem Registration
// ============================================================================
void RussianDollProblem::setPVProblem(CoupledProblem* pvProblem) {
    M_pvProblem = pvProblem;
    std::cout << "[RussianDoll] PV problem added" << std::endl;
}

void RussianDollProblem::setPLCProblem(PLCProblem* plcProblem) {
    M_plcProblem = plcProblem;
    
    // Set the outer wall region
    plcProblem->setOuterWallRegion(M_outerWallRegion);
    
    std::cout << "[RussianDoll] PLC problem added" << std::endl;
}

// ============================================================================
// Initialize
// ============================================================================
void RussianDollProblem::initialize() {
    if (M_initialized) {
        std::cout << "[RussianDoll] Already initialized" << std::endl;
        return;
    }
    
    if (!M_pvProblem || !M_plcProblem) {
        std::cerr << "[RussianDoll] Error: Problems not registered!" << std::endl;
        return;
    }
    
    // Initialize interpolation manager
    M_interpManager->initialize();
    
    // Get z-range from PV mesh
    getfem::mesh* pvMesh = M_bulkPV->getMesh();
    M_z_min = std::numeric_limits<scalar_type>::max();
    M_z_max = std::numeric_limits<scalar_type>::lowest();
    
    for (dal::bv_visitor i(pvMesh->points_index()); !i.finished(); ++i) {
        scalar_type z = pvMesh->points()[i][2];
        if (z < M_z_min) M_z_min = z;
        if (z > M_z_max) M_z_max = z;
    }
    
    std::cout << "[RussianDoll] Z-range from PV mesh: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    
    // Initialize coefficient vectors
    M_pressureCoefficients.resize(M_polynomialOrder + 1, 0.0);
    M_displacementXCoefficients.resize(M_polynomialOrder + 1, 0.0);
    M_displacementYCoefficients.resize(M_polynomialOrder + 1, 0.0);
    M_displacementZCoefficients.resize(M_polynomialOrder + 1, 0.0);
    
    M_initialized = true;
    std::cout << "[RussianDoll] Initialization complete" << std::endl;
}

// ============================================================================
// Core Coupling: Interpolate PV to PLC
// ============================================================================
void RussianDollProblem::interpolatePVtoPLC() {
    std::cout << "\n--- RussianDoll: Interpolating PV -> PLC ---" << std::endl;
    
    if (!M_initialized) {
        std::cerr << "[RussianDoll] Error: Not initialized!" << std::endl;
        return;
    }
    
    if (!M_pvProblem || !M_plcProblem) {
        std::cerr << "[RussianDoll] Error: Problems not set!" << std::endl;
        return;
    }
    
    // ========================================================================
    // Step 1: Get PV solutions
    // ========================================================================
    scalarVectorPtr_Type pvPressure = M_pvProblem->getPressure();
    scalarVectorPtr_Type pvDisplacement = M_pvProblem->getDisplacement();
    
    if (!pvPressure || pvPressure->empty()) {
        std::cerr << "[RussianDoll] Error: PV pressure solution not available!" << std::endl;
        return;
    }
    #ifdef VERBOSE
    std::cout << "  PV pressure DOFs: " << pvPressure->size() << std::endl;
    std::cout << "  PV pressure norm: " << gmm::vect_norm2(*pvPressure) << std::endl;
 

    if (pvDisplacement && !pvDisplacement->empty()) {
        std::cout << "  PV displacement DOFs: " << pvDisplacement->size() << std::endl;
        std::cout << "  PV displacement norm: " << gmm::vect_norm2(*pvDisplacement) << std::endl;
    }
   #endif    
    // ========================================================================
    // Step 2: Use InterpolationManager to extract and fit
    // ========================================================================
    const getfem::mesh_fem& mf_pv_pressure = M_pvProblem->getMfPressure();
    
    // Interpolate pressure
    PolynomialFit pressureFit = M_interpManager->interpolatePressure(
        pvPressure, mf_pv_pressure);
    
    // Store coefficients
    M_pressureCoefficients = pressureFit.coefficients;
    M_z_min = pressureFit.z_min;
    M_z_max = pressureFit.z_max;
    
    std::cout << "  Pressure polynomial fit:" << std::endl;
    std::cout << "    Z-range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    std::cout << "    R² = " << pressureFit.r_squared << std::endl;
    std::cout << "    Coefficients: [";
    for (size_t i = 0; i < M_pressureCoefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << std::scientific << std::setprecision(4) << M_pressureCoefficients[i];
    }
    std::cout << "]" << std::endl;
    
    // Interpolate displacement if available
    if (pvDisplacement && !pvDisplacement->empty()) {
        const getfem::mesh_fem& mf_pv_disp = M_pvProblem->getMfDisplacement();
        
        PolynomialFit fitX, fitY, fitZ;
        M_interpManager->interpolateDisplacement(
            pvDisplacement, mf_pv_disp, fitX, fitY, fitZ);
        
        M_displacementXCoefficients = fitX.coefficients;
        M_displacementYCoefficients = fitY.coefficients;
        M_displacementZCoefficients = fitZ.coefficients;
        
        std::cout << "  Displacement polynomial fit completed" << std::endl;
    }
    
    // ========================================================================
    // Step 3: Update PLC boundary conditions
    // ========================================================================
    
    // Set pressure BC (Neumann for Darcy in mixed form)
    // This sets the callback in BC class which is used by naturalRHS()
    M_plcProblem->setOuterWallPressureBC(M_pressureCoefficients, M_z_min, M_z_max);
    
    // Set displacement BC (Dirichlet for elasticity)
    // This sets the callback in BC class which is used by BCDiriVec()
    if (!M_displacementXCoefficients.empty()) {
        M_plcProblem->setOuterWallDisplacementBC(
            M_displacementXCoefficients,
            M_displacementYCoefficients,
            M_displacementZCoefficients,
            M_z_min, M_z_max);
    }
    
    std::cout << "--- Interpolation complete ---\n" << std::endl;
}

// ============================================================================
// Update Solutions
// ============================================================================
void RussianDollProblem::updateSolutions() {
    if (M_pvProblem) {
        M_pvProblem->updateSol();
    }
    if (M_plcProblem) {
        M_plcProblem->updateSol();
    }
}

// ============================================================================
// Export VTK
// ============================================================================
void RussianDollProblem::exportVtk(const std::string& folder, int frame) {
    std::cout << "[RussianDoll] Exporting solutions (frame " << frame << ")..." << std::endl;
    
    if (M_pvProblem) {
        M_pvProblem->exportVtk(folder + "/pv", "all", frame);
    }
    
    if (M_plcProblem) {
        M_plcProblem->exportVtk(folder + "/plc", "all", frame);
    }
}

// ============================================================================
// Export Interpolation Data
// ============================================================================
void RussianDollProblem::exportInterpolationData(const std::string& folder, int frame) {
    std::string filename = folder + "/interpolation_data_" + std::to_string(frame) + ".txt";
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[RussianDoll] Error: Cannot open file " << filename << std::endl;
        return;
    }
    
    file << "# RussianDoll Interpolation Data - Frame " << frame << std::endl;
    file << "# Time: " << M_time->time() << std::endl;
    file << std::endl;
    
    file << "# Z-range" << std::endl;
    file << "z_min = " << M_z_min << std::endl;
    file << "z_max = " << M_z_max << std::endl;
    file << std::endl;
    
    file << "# Pressure polynomial coefficients (order " << (M_pressureCoefficients.size()-1) << ")" << std::endl;
    file << "# p(t) = c0 + c1*t + c2*t^2 + ... where t = (z - z_min)/(z_max - z_min)" << std::endl;
    for (size_t i = 0; i < M_pressureCoefficients.size(); ++i) {
        file << "c" << i << " = " << std::scientific << std::setprecision(10) 
             << M_pressureCoefficients[i] << std::endl;
    }
    file << std::endl;
    
    if (!M_displacementXCoefficients.empty()) {
        file << "# Displacement X coefficients" << std::endl;
        for (size_t i = 0; i < M_displacementXCoefficients.size(); ++i) {
            file << "ux_c" << i << " = " << std::scientific << std::setprecision(10)
                 << M_displacementXCoefficients[i] << std::endl;
        }
        file << std::endl;
        
        file << "# Displacement Y coefficients" << std::endl;
        for (size_t i = 0; i < M_displacementYCoefficients.size(); ++i) {
            file << "uy_c" << i << " = " << std::scientific << std::setprecision(10)
                 << M_displacementYCoefficients[i] << std::endl;
        }
        file << std::endl;
        
        file << "# Displacement Z coefficients" << std::endl;
        for (size_t i = 0; i < M_displacementZCoefficients.size(); ++i) {
            file << "uz_c" << i << " = " << std::scientific << std::setprecision(10)
                 << M_displacementZCoefficients[i] << std::endl;
        }
    }
    file << std::endl;
    
    // Also write sampled values for plotting
    file << "# Sampled values along z-axis" << std::endl;
    file << "# z, p_v(z), ux(z), uy(z), uz(z)" << std::endl;
    
    int nSamples = 50;
    for (int i = 0; i <= nSamples; ++i) {
        scalar_type t = static_cast<scalar_type>(i) / nSamples;
        scalar_type z = M_z_min + t * (M_z_max - M_z_min);
        
        scalar_type p = evaluatePolynomial(M_pressureCoefficients, z);
        scalar_type ux = evaluatePolynomial(M_displacementXCoefficients, z);
        scalar_type uy = evaluatePolynomial(M_displacementYCoefficients, z);
        scalar_type uz = evaluatePolynomial(M_displacementZCoefficients, z);
        
        file << std::fixed << std::setprecision(6) << z << ", "
             << std::scientific << std::setprecision(6) << p << ", "
             << ux << ", " << uy << ", " << uz << std::endl;
    }
    
    file.close();
    std::cout << "[RussianDoll] Interpolation data exported to " << filename << std::endl;
}

// ============================================================================
// Compute Errors
// ============================================================================
std::vector<scalar_type> RussianDollProblem::computeErrors(scalar_type time) {
    std::vector<scalar_type> errors(4, -1.0);  // [pv_err, uv_err, pl_err, ul_err]
    
    if (M_pvProblem) {
        bgeot::base_node pvErrors = M_pvProblem->computeError(time);
        errors[0] = pvErrors[0];  // Pressure error
        errors[1] = pvErrors[1];  // Displacement error
    }
    
    if (M_plcProblem) {
        bgeot::base_node plcErrors = M_plcProblem->computeError(time);
        errors[2] = plcErrors[0];  // Pressure error
        errors[3] = plcErrors[1];  // Displacement error
    }
    
    return errors;
}

// ============================================================================
// Print Coupling Info
// ============================================================================
void RussianDollProblem::printCouplingInfo() const {
    std::cout << "\n=====================================================" << std::endl;
    std::cout << " RussianDoll Coupling Information" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Initialized: " << (M_initialized ? "YES" : "NO") << std::endl;
    std::cout << "  Polynomial order: " << M_polynomialOrder << std::endl;
    std::cout << "  Z-range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    std::cout << "  Outer wall region: " << M_outerWallRegion << std::endl;
    
    if (!M_pressureCoefficients.empty()) {
        std::cout << "  Pressure coefficients: [";
        for (size_t i = 0; i < M_pressureCoefficients.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << std::scientific << std::setprecision(3) << M_pressureCoefficients[i];
        }
        std::cout << "]" << std::endl;
    }
    
    if (M_plcProblem) {
        M_plcProblem->printCouplingInfo();
    }
    
    std::cout << "=====================================================" << std::endl;
}

// ============================================================================
// Get Displacement Coefficients
// ============================================================================
const std::vector<scalar_type>& RussianDollProblem::getDisplacementCoefficients(int component) const {
    switch (component) {
        case 0: return M_displacementXCoefficients;
        case 1: return M_displacementYCoefficients;
        case 2: return M_displacementZCoefficients;
        default:
            std::cerr << "[RussianDoll] Warning: Invalid component " << component << std::endl;
            return M_displacementXCoefficients;
    }
}

// ============================================================================
// Evaluate Polynomial
// ============================================================================
scalar_type RussianDollProblem::evaluatePolynomial(const std::vector<scalar_type>& coeffs, 
                                                    scalar_type z) const {
    if (coeffs.empty()) return 0.0;
    
    // Normalize z to [0, 1]
    scalar_type t = 0.0;
    if (std::abs(M_z_max - M_z_min) > 1e-15) {
        t = (z - M_z_min) / (M_z_max - M_z_min);
    }
    t = std::max(0.0, std::min(1.0, t));
    
    // Evaluate: p(t) = c0 + c1*t + c2*t^2 + ...
    scalar_type result = 0.0;
    scalar_type t_power = 1.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * t_power;
        t_power *= t;
    }
    
    return result;
}