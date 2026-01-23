// ============================================================================
// RussianDollProblem.cc - Implementation of dual-porosity coupling
// ============================================================================
// This version implements one-way coupling from PV to PLC with polynomial
// BC coefficients extracted via line interpolation.
// ============================================================================

#include "../include/RussianDollProblem.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

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
      M_gamma(0.0),  // One-way coupling: no feedback from PLC to PV
      M_z_min(0.0),
      M_z_max(1.0),
      M_section("russian_doll/"),
      M_initialized(false),
      M_polynomialOrder(3)
{
    // Read coupling parameters from data file
    M_gamma = dataFile((M_section + "leakage_coefficient").c_str(), 0.0);
    M_polynomialOrder = dataFile((M_section + "polynomial_order").c_str(), 3);
    
    // Read coupling approach
    std::string approach_str = dataFile((M_section + "interpolation/approach").c_str(), "line");
    if (approach_str == "mesh" || approach_str == "MESH") {
        M_couplingApproach = CouplingApproach::MESH_INTERPOLATION;
    } else {
        M_couplingApproach = CouplingApproach::LINE_INTERPOLATION;
    }
    
    // Create interpolation manager
    M_interpManager = std::make_unique<InterpolationManager>(dataFile, bulkPV, bulkPLC);
    
    // Initialize coefficient vector
    M_pvBCCoefficients.resize(M_polynomialOrder + 1, 0.0);
    
    std::cout << "=====================================================" << std::endl;
    std::cout << " RussianDollProblem Configuration (One-Way Coupling)" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << " Leakage coefficient (gamma): " << M_gamma << std::endl;
    std::cout << " Polynomial order: " << M_polynomialOrder << std::endl;
    std::cout << " Coupling approach: " 
              << (M_couplingApproach == CouplingApproach::LINE_INTERPOLATION 
                  ? "LINE_INTERPOLATION" : "MESH_INTERPOLATION") << std::endl;
    std::cout << "=====================================================" << std::endl;
}

// ============================================================================
// Setup Methods
// ============================================================================
void RussianDollProblem::setPVProblem(CoupledProblem* pvProblem) {
    M_pvProblem = pvProblem;
    std::cout << "[RussianDoll] PV problem registered" << std::endl;
}

void RussianDollProblem::setPLCProblem(PLCProblem* plcProblem) {
    M_plcProblem = plcProblem;
    
    // Set leakage coefficient if doing two-way coupling (typically 0 for one-way)
    M_plcProblem->setLeakageCoefficient(M_gamma);
    
    std::cout << "[RussianDoll] PLC problem registered" << std::endl;
}

void RussianDollProblem::setupLineProfile() {
    // Setup a vertical line profile through the domain center
    // This line is used to extract PV pressure vs z
    
    // Get domain bounds from PV mesh
    const getfem::mesh& mesh = *(M_bulkPV->getMesh());
    
    // Find bounding box
    scalar_type x_min = std::numeric_limits<scalar_type>::max();
    scalar_type x_max = std::numeric_limits<scalar_type>::lowest();
    scalar_type y_min = std::numeric_limits<scalar_type>::max();
    scalar_type y_max = std::numeric_limits<scalar_type>::lowest();
    M_z_min = std::numeric_limits<scalar_type>::max();
    M_z_max = std::numeric_limits<scalar_type>::lowest();
    
    for (dal::bv_visitor i(mesh.points().index()); !i.finished(); ++i) {
        bgeot::base_node pt = mesh.points()[i];
        if (pt.size() >= 3) {
            x_min = std::min(x_min, pt[0]);
            x_max = std::max(x_max, pt[0]);
            y_min = std::min(y_min, pt[1]);
            y_max = std::max(y_max, pt[1]);
            M_z_min = std::min(M_z_min, pt[2]);
            M_z_max = std::max(M_z_max, pt[2]);
        }
    }
    
    // Create vertical line through domain center
    scalar_type x_center = (x_min + x_max) / 2.0;
    scalar_type y_center = (y_min + y_max) / 2.0;
    
    LineProfile profile;
    profile.start_point.resize(3);
    profile.start_point[0] = x_center;
    profile.start_point[1] = y_center;
    profile.start_point[2] = M_z_min;
    
    profile.end_point.resize(3);
    profile.end_point[0] = x_center;
    profile.end_point[1] = y_center;
    profile.end_point[2] = M_z_max;
    
    profile.num_samples = 100;  // Adequate for polynomial fitting
    profile.polynomial_order = M_polynomialOrder;
    profile.name = "pv_vertical_line";
    
    // Add to interpolation manager
    M_interpManager->addLineProfile(profile);
    
    std::cout << "[RussianDoll] Line profile created:" << std::endl;
    std::cout << "  Start: (" << profile.start_point[0] << ", " 
              << profile.start_point[1] << ", " << profile.start_point[2] << ")" << std::endl;
    std::cout << "  End: (" << profile.end_point[0] << ", " 
              << profile.end_point[1] << ", " << profile.end_point[2] << ")" << std::endl;
    std::cout << "  Z range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
}

void RussianDollProblem::initialize() {
    if (M_pvProblem == nullptr || M_plcProblem == nullptr) {
        std::cerr << "[RussianDoll] Error: PV or PLC problem not set!" << std::endl;
        return;
    }
    
    // Setup line profile for interpolation
    if (M_couplingApproach == CouplingApproach::LINE_INTERPOLATION) {
        setupLineProfile();
    }
    
    // Initialize interpolation manager
    M_interpManager->initialize();
    
    // Setup BC callback in PLC problem
    // Capture 'this' to access polynomial evaluation
    M_plcProblem->setOuterWallBCCallback(
        [this](scalar_type z) -> scalar_type {
            return this->evaluatePVPressure(z);
        }
    );
    
    // Set z-range in PLC for proper normalization
    // (Will be updated after first interpolation)
    
    M_initialized = true;
    
    std::cout << "[RussianDoll] Initialization complete" << std::endl;
    std::cout << "  PV pressure DOFs: " << M_pvProblem->getNbPressureDOF() << std::endl;
    std::cout << "  PLC pressure DOFs: " << M_plcProblem->getNbPressureDOF() << std::endl;
}

// ============================================================================
// Solution Methods
// ============================================================================
void RussianDollProblem::solveTimeStep() {
    if (!M_initialized) {
        std::cerr << "[RussianDoll] Error: Not initialized!" << std::endl;
        return;
    }
    
    std::cout << "\n=== RussianDoll Time Step ===" << std::endl;
    
    // Step 1: Solve PV problem (standalone, one-way coupling)
    solvePV();
    
    // Step 2: Extract PV pressure and fit polynomial
    interpolatePVtoPLC();
    
    // Step 3: Update PLC boundary coefficients
    updatePLCBoundaryCoefficients();
    
    // Step 4: Solve PLC problem with polynomial BC
    solvePLC();
    
    std::cout << "=== Time Step Complete ===" << std::endl;
}

void RussianDollProblem::solvePV() {
    std::cout << "[RussianDoll] Solving PV problem..." << std::endl;
    
    // PV is solved standalone - no coupling from PLC (one-way)
    // The assembly and solve are handled by the main time loop
    // This method is called after PV has been assembled
    
    // Note: The actual solve happens in the main time loop
    // Here we just mark that PV solution is ready for interpolation
}

void RussianDollProblem::solvePLC() {
    std::cout << "[RussianDoll] Solving PLC problem..." << std::endl;
    
    // PLC uses the polynomial BC set via updatePLCBoundaryCoefficients()
    // The actual solve happens in the main time loop
}

void RussianDollProblem::interpolatePVtoPLC() {
    std::cout << "[RussianDoll] Interpolating PV pressure to polynomial..." << std::endl;
    
    // Get PV pressure solution
    scalarVectorPtr_Type pvPressure = M_pvProblem->getPressure();
    
    if (!pvPressure || pvPressure->size() == 0) {
        std::cerr << "[RussianDoll] Warning: PV pressure not available!" << std::endl;
        return;
    }
    
    if (M_couplingApproach == CouplingApproach::LINE_INTERPOLATION) {
        // Use line interpolation approach
        const getfem::mesh_fem& mf_pv = M_pvProblem->getMfPressure();
        
        // Get line profiles
        const std::vector<LineProfile>& profiles = M_interpManager->getLineProfiles();
        
        if (profiles.empty()) {
            std::cerr << "[RussianDoll] Error: No line profiles defined!" << std::endl;
            return;
        }
        
        // Use first line profile (vertical line)
        const LineProfile& profile = profiles[0];
        
        // Extract values along line
        std::vector<scalar_type> arc_coords, values;
        M_interpManager->extractAlongLine(pvPressure, mf_pv, profile, arc_coords, values);
        
        // Fit polynomial
        PolynomialFit fit = M_interpManager->fitPolynomial(arc_coords, values, profile.polynomial_order);
        
        // Store coefficients
        M_pvBCCoefficients = fit.coefficients;
        
        std::cout << "[RussianDoll] Polynomial fit complete:" << std::endl;
        std::cout << "  Order: " << profile.polynomial_order << std::endl;
        std::cout << "  R² = " << fit.r_squared << std::endl;
        std::cout << "  Coefficients: [";
        for (size_t i = 0; i < M_pvBCCoefficients.size(); ++i) {
            std::cout << M_pvBCCoefficients[i];
            if (i < M_pvBCCoefficients.size() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
        
    } else {
        // Mesh interpolation approach - not typical for BC transfer
        // Would need different handling (direct interpolation to all PLC DOFs)
        std::cerr << "[RussianDoll] Warning: MESH_INTERPOLATION not ideal for BC transfer" << std::endl;
    }
}

void RussianDollProblem::updatePLCBoundaryCoefficients() {
    std::cout << "[RussianDoll] Updating PLC boundary coefficients..." << std::endl;
    
    // Pass polynomial coefficients to PLC problem
    M_plcProblem->setOuterWallBCCoefficients(M_pvBCCoefficients, M_z_min, M_z_max);
    
    // The BC callback was already set in initialize(), so PLC can evaluate p_v(z)
}

void RussianDollProblem::updateSolutions() {
    std::cout << "[RussianDoll] Updating solutions for next time step..." << std::endl;
    
    // Update is handled by individual problems in main time loop
    M_pvProblem->updateSol();
    M_plcProblem->updateSol();
}

// ============================================================================
// Coefficient Access
// ============================================================================
scalar_type RussianDollProblem::evaluatePVPressure(scalar_type z) const {
    if (M_pvBCCoefficients.empty()) {
        return 0.0;
    }
    
    // Normalize z to [0, 1]
    scalar_type t = 0.0;
    if (std::abs(M_z_max - M_z_min) > 1e-15) {
        t = (z - M_z_min) / (M_z_max - M_z_min);
    }
    
    // Clamp to [0, 1]
    t = std::max(0.0, std::min(1.0, t));
    
    // Evaluate polynomial: p(t) = c0 + c1*t + c2*t^2 + ...
    scalar_type result = 0.0;
    scalar_type t_power = 1.0;
    
    for (size_type i = 0; i < M_pvBCCoefficients.size(); ++i) {
        result += M_pvBCCoefficients[i] * t_power;
        t_power *= t;
    }
    
    return result;
}

void RussianDollProblem::getZRange(scalar_type& z_min, scalar_type& z_max) const {
    z_min = M_z_min;
    z_max = M_z_max;
}

// ============================================================================
// Output Methods
// ============================================================================
void RussianDollProblem::exportVtk(const std::string& folder, int frame) {
    std::cout << "[RussianDoll] Exporting VTK (frame " << frame << ")..." << std::endl;
    
    // Export PV solution
    M_pvProblem->exportVtk(folder + "/pv", "all", frame);
    
    // Export PLC solution
    M_plcProblem->exportVtk(folder + "/plc", "all", frame);
}

void RussianDollProblem::exportInterpolationData(const std::string& folder, int frame) {
    // Export polynomial coefficients and fit data for debugging
    std::string filename = folder + "/interpolation_data_" + std::to_string(frame) + ".txt";
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[RussianDoll] Could not open " << filename << std::endl;
        return;
    }
    
    file << "# RussianDoll Interpolation Data - Frame " << frame << std::endl;
    file << "# Z range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    file << "# Polynomial order: " << M_polynomialOrder << std::endl;
    file << "# Coefficients (c0, c1, c2, ...):" << std::endl;
    
    for (size_type i = 0; i < M_pvBCCoefficients.size(); ++i) {
        file << "c" << i << " = " << M_pvBCCoefficients[i] << std::endl;
    }
    
    file << std::endl;
    file << "# Polynomial values at sampled z:" << std::endl;
    file << "# z\tp_v(z)" << std::endl;
    
    size_type num_samples = 50;
    for (size_type i = 0; i <= num_samples; ++i) {
        scalar_type t = static_cast<scalar_type>(i) / num_samples;
        scalar_type z = M_z_min + t * (M_z_max - M_z_min);
        scalar_type p = evaluatePVPressure(z);
        file << z << "\t" << p << std::endl;
    }
    
    file.close();
    std::cout << "[RussianDoll] Interpolation data exported to " << filename << std::endl;
}

std::vector<scalar_type> RussianDollProblem::computeErrors(scalar_type time) {
    std::vector<scalar_type> errors(4, 0.0);
    
    // Compute PV errors
    bgeot::base_node pvErrors = M_pvProblem->computeError(time);
    errors[0] = pvErrors[0];  // PV pressure error
    errors[1] = pvErrors[1];  // PV displacement error
    
    // Compute PLC errors
    bgeot::base_node plcErrors = M_plcProblem->computeError(time);
    errors[2] = plcErrors[0];  // PLC pressure error
    errors[3] = plcErrors[1];  // PLC displacement error
    
    std::cout << "[RussianDoll] Errors at t=" << time << ":" << std::endl;
    std::cout << "  PV:  p_error=" << errors[0] << ", u_error=" << errors[1] << std::endl;
    std::cout << "  PLC: p_error=" << errors[2] << ", u_error=" << errors[3] << std::endl;
    
    return errors;
}