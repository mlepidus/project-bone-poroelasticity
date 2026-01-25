// ============================================================================
// RussianDollProblem.cc - Implementation of dual-porosity coupling
// ============================================================================
// Complete workflow:
// 1. Solve PV problem
// 2. Extract PV pressure/displacement along line
// 3. Fit polynomials and store coefficients
// 4. Update PLC boundary conditions
// 5. Compute coupling RHS: +gamma * M * p_v
// 6. Solve PLC problem
// ============================================================================

#include "../include/RussianDollProblem.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
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
      M_z_min(0.0),
      M_z_max(1.0),
      M_gamma(0.0),
      M_section("russian_doll/"),
      M_initialized(false),
      M_polynomialOrder(3),
      M_pvPressureOnPLC(nullptr)
{
    // Read coupling parameters from data file
    M_gamma = dataFile((M_section + "leakage_coefficient").c_str(), 
                       dataFile("bulkDataPLC/darcy/leakage", 0.0));
    M_polynomialOrder = dataFile((M_section + "polynomial_order").c_str(), 
                                  dataFile((M_section + "interpolation/line0/polynomial_order").c_str(), 3));
    
    // Read coupling approach
    std::string approach_str = dataFile((M_section + "interpolation/approach").c_str(), "line");
    if (approach_str == "mesh" || approach_str == "MESH") {
        M_couplingApproach = CouplingApproach::MESH_INTERPOLATION;
    } else {
        M_couplingApproach = CouplingApproach::LINE_INTERPOLATION;
    }
    
    // Create interpolation manager
    M_interpManager = std::make_unique<InterpolationManager>(dataFile, bulkPV, bulkPLC);
    
    // Initialize coefficient vectors
    M_pressureCoefficients.resize(M_polynomialOrder + 1, 0.0);
    M_displacementXCoefficients.resize(M_polynomialOrder + 1, 0.0);
    M_displacementYCoefficients.resize(M_polynomialOrder + 1, 0.0);
    M_displacementZCoefficients.resize(M_polynomialOrder + 1, 0.0);
    
    std::cout << "=====================================================" << std::endl;
    std::cout << " RussianDollProblem Configuration" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << " Leakage coefficient (gamma): " << M_gamma << std::endl;
    std::cout << " Polynomial order: " << M_polynomialOrder << std::endl;
    std::cout << " Interpolation approach: " 
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
    
    // Set leakage coefficient in PLC problem
    M_plcProblem->setLeakageCoefficient(M_gamma);
    
    std::cout << "[RussianDoll] PLC problem registered with gamma = " << M_gamma << std::endl;
}

void RussianDollProblem::setupLineProfile() {
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
    
    profile.num_samples = 100;
    profile.polynomial_order = M_polynomialOrder;
    profile.name = "pv_vertical_line";
    
    // Add to interpolation manager if not already configured
    if (M_interpManager->getLineProfiles().empty()) {
        M_interpManager->addLineProfile(profile);
    }
    
    std::cout << "[RussianDoll] Line profile configured:" << std::endl;
    std::cout << "  Center: (" << x_center << ", " << y_center << ")" << std::endl;
    std::cout << "  Z range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    std::cout << "  Polynomial order: " << M_polynomialOrder << std::endl;
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
    M_plcProblem->setOuterWallBCCallback(
        [this](scalar_type z) -> scalar_type {
            return this->evaluatePVPressure(z);
        }
    );
    
    // Allocate PV pressure vector for PLC DOFs
    size_type nbPressureDOF = M_plcProblem->getNbPressureDOF();
    M_pvPressureOnPLC.reset(new scalarVector_Type(nbPressureDOF, 0.0));
    
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
    std::cout << "Step 1: Solving PV problem..." << std::endl;
    solvePV();
    
    // Step 2: Extract PV pressure and fit polynomial
    std::cout << "Step 2: Interpolating PV -> PLC..." << std::endl;
    interpolatePVtoPLC();
    
    // Step 3: Update PLC boundary coefficients and coupling source
    std::cout << "Step 3: Updating PLC boundary and coupling..." << std::endl;
    updatePLCBoundaryCoefficients();
    
    // Step 4: Solve PLC problem with polynomial BC
    std::cout << "Step 4: Solving PLC problem..." << std::endl;
    solvePLC();
    
    std::cout << "=== Time Step Complete ===" << std::endl;
}

void RussianDollProblem::solvePV() {
    // PV is solved standalone - no coupling from PLC (one-way)
    // The actual assembly and solve are handled by the main time loop
    std::cout << "[RussianDoll] PV solve marked complete" << std::endl;
}

void RussianDollProblem::solvePLC() {
    // PLC uses the polynomial BC set via updatePLCBoundaryCoefficients()
    // The actual assembly and solve are handled by the main time loop
    std::cout << "[RussianDoll] PLC solve marked complete" << std::endl;
}

void RussianDollProblem::interpolatePVtoPLC() {
    std::cout << "[RussianDoll] Starting interpolation..." << std::endl;
    
    // Get PV pressure solution
    scalarVectorPtr_Type pvPressure = M_pvProblem->getPressure();
    
    if (!pvPressure || pvPressure->size() == 0) {
        std::cerr << "[RussianDoll] Warning: PV pressure not available!" << std::endl;
        return;
    }
    
    std::cout << "  PV pressure solution: " << pvPressure->size() << " DOFs, "
              << "norm = " << gmm::vect_norm2(*pvPressure) << std::endl;
    
    if (M_couplingApproach == CouplingApproach::LINE_INTERPOLATION) {
        // Use line interpolation approach
        const getfem::mesh_fem& mf_pv = M_pvProblem->getMfPressure();
        
        // Interpolate pressure
        PolynomialFit pressureFit = M_interpManager->interpolatePressure(pvPressure, mf_pv);
        
        // Store pressure coefficients
        M_pressureCoefficients = pressureFit.coefficients;
        M_z_min = pressureFit.z_min;
        M_z_max = pressureFit.z_max;
        
        // Optionally interpolate displacement
        scalarVectorPtr_Type pvDisplacement = M_pvProblem->getDisplacement();
        if (pvDisplacement && pvDisplacement->size() > 0) {
            const getfem::mesh_fem& mf_disp = M_pvProblem->getMfDisplacement();
            
            PolynomialFit fit_x, fit_y, fit_z;
            M_interpManager->interpolateDisplacement(pvDisplacement, mf_disp,
                                                     fit_x, fit_y, fit_z);
            
            M_displacementXCoefficients = fit_x.coefficients;
            M_displacementYCoefficients = fit_y.coefficients;
            M_displacementZCoefficients = fit_z.coefficients;
        }
        
        // Print coefficient information
        printCoefficients();
        
    } else {
        // Mesh interpolation approach
        std::cerr << "[RussianDoll] Warning: MESH_INTERPOLATION requires different handling" << std::endl;
    }
}

void RussianDollProblem::updatePLCBoundaryCoefficients() {
    std::cout << "[RussianDoll] Updating PLC boundary coefficients..." << std::endl;
    
    // Pass polynomial coefficients to PLC problem
    M_plcProblem->setOuterWallBCCoefficients(M_pressureCoefficients, M_z_min, M_z_max);
    
    M_plcProblem->setOuterWallDisplacementBCCoefficients(M_displacementXCoefficients,
                                                        M_displacementYCoefficients,
                                                        M_displacementZCoefficients,
                                                        M_z_min, M_z_max);
    // Compute PV pressure at PLC DOF locations for coupling RHS
    computePVPressureOnPLCDOFs();
    
    // Set the coupling source in PLC problem
    M_plcProblem->setCouplingSource(M_pvPressureOnPLC);
    
    std::cout << "[RussianDoll] PLC coupling updated" << std::endl;
}

void RussianDollProblem::computePVPressureOnPLCDOFs() {
    if (!M_plcProblem || !M_plcProblem->getDarcyPB()) {
        std::cerr << "[RussianDoll] Error: PLC Darcy problem not available!" << std::endl;
        return;
    }
    
    // Get PLC pressure mesh_fem
    const getfem::mesh_fem& mf_plc_pressure = 
        *(M_plcProblem->getDarcyPB()->getFEM("Pressure")->getFEM());
    
    size_type nbPressureDOF = mf_plc_pressure.nb_dof();
    
    // Resize if needed
    if (M_pvPressureOnPLC->size() != nbPressureDOF) {
        M_pvPressureOnPLC.reset(new scalarVector_Type(nbPressureDOF, 0.0));
    }
    
    // Evaluate polynomial at each PLC pressure DOF location
    for (size_type i = 0; i < nbPressureDOF; ++i) {
        bgeot::base_node pt = mf_plc_pressure.point_of_basic_dof(i);
        scalar_type z = pt[2];
        (*M_pvPressureOnPLC)[i] = evaluatePVPressure(z);
    }
    
    std::cout << "[RussianDoll] Computed PV pressure at " << nbPressureDOF 
              << " PLC DOFs (norm = " << gmm::vect_norm2(*M_pvPressureOnPLC) << ")" << std::endl;
}

void RussianDollProblem::updateSolutions() {
    std::cout << "[RussianDoll] Updating solutions for next time step..." << std::endl;
    M_pvProblem->updateSol();
    M_plcProblem->updateSol();
}

// ============================================================================
// Coefficient Access
// ============================================================================
const std::vector<scalar_type>& RussianDollProblem::getDisplacementCoefficients(int component) const {
    switch (component) {
        case 0: return M_displacementXCoefficients;
        case 1: return M_displacementYCoefficients;
        case 2: return M_displacementZCoefficients;
        default: 
            static std::vector<scalar_type> empty;
            return empty;
    }
}

scalar_type RussianDollProblem::evaluatePVPressure(scalar_type z) const {
    if (M_pressureCoefficients.empty()) {
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
    
    for (size_type i = 0; i < M_pressureCoefficients.size(); ++i) {
        result += M_pressureCoefficients[i] * t_power;
        t_power *= t;
    }
    
    return result;
}

void RussianDollProblem::evaluatePVDisplacement(scalar_type z, 
                                                 scalar_type& ux, 
                                                 scalar_type& uy, 
                                                 scalar_type& uz) const {
    // Normalize z to [0, 1]
    scalar_type t = 0.0;
    if (std::abs(M_z_max - M_z_min) > 1e-15) {
        t = (z - M_z_min) / (M_z_max - M_z_min);
    }
    t = std::max(0.0, std::min(1.0, t));
    
    // Evaluate each component
    ux = 0.0; uy = 0.0; uz = 0.0;
    scalar_type t_power = 1.0;
    
    for (size_type i = 0; i < M_displacementXCoefficients.size(); ++i) {
        ux += M_displacementXCoefficients[i] * t_power;
        uy += M_displacementYCoefficients[i] * t_power;
        uz += M_displacementZCoefficients[i] * t_power;
        t_power *= t;
    }
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
    
    file << std::setprecision(12);
    file << "# RussianDoll Interpolation Data - Frame " << frame << std::endl;
    file << "# =================================================" << std::endl;
    file << "# Z range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    file << "# Polynomial order: " << M_polynomialOrder << std::endl;
    file << "# Leakage coefficient (gamma): " << M_gamma << std::endl;
    file << "#" << std::endl;
    
    // Pressure coefficients
    file << "# Pressure polynomial coefficients p_v(t) = sum_i c_i * t^i" << std::endl;
    file << "# where t = (z - z_min) / (z_max - z_min)" << std::endl;
    file << "[pressure_coefficients]" << std::endl;
    for (size_type i = 0; i < M_pressureCoefficients.size(); ++i) {
        file << "c" << i << " = " << M_pressureCoefficients[i] << std::endl;
    }
    file << std::endl;
    
    // Displacement coefficients
    file << "[displacement_x_coefficients]" << std::endl;
    for (size_type i = 0; i < M_displacementXCoefficients.size(); ++i) {
        file << "c" << i << " = " << M_displacementXCoefficients[i] << std::endl;
    }
    file << std::endl;
    
    file << "[displacement_y_coefficients]" << std::endl;
    for (size_type i = 0; i < M_displacementYCoefficients.size(); ++i) {
        file << "c" << i << " = " << M_displacementYCoefficients[i] << std::endl;
    }
    file << std::endl;
    
    file << "[displacement_z_coefficients]" << std::endl;
    for (size_type i = 0; i < M_displacementZCoefficients.size(); ++i) {
        file << "c" << i << " = " << M_displacementZCoefficients[i] << std::endl;
    }
    file << std::endl;
    
    // Sampled values
    file << "# Polynomial values at sampled z:" << std::endl;
    file << "# z\tt\tp_v(z)\tu_x(z)\tu_y(z)\tu_z(z)" << std::endl;
    file << "[sampled_values]" << std::endl;
    
    size_type num_samples = 50;
    for (size_type i = 0; i <= num_samples; ++i) {
        scalar_type t_val = static_cast<scalar_type>(i) / num_samples;
        scalar_type z = M_z_min + t_val * (M_z_max - M_z_min);
        scalar_type p = evaluatePVPressure(z);
        scalar_type ux, uy, uz;
        evaluatePVDisplacement(z, ux, uy, uz);
        
        file << z << "\t" << t_val << "\t" << p << "\t" 
             << ux << "\t" << uy << "\t" << uz << std::endl;
    }
    
    file.close();
    std::cout << "[RussianDoll] Interpolation data exported to " << filename << std::endl;
}

void RussianDollProblem::exportCoefficients(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[RussianDoll] Could not open " << filename << std::endl;
        return;
    }
    
    file << std::setprecision(15);
    
    // Export in a format that can be easily read back
    file << "# Polynomial coefficients for Russian Doll coupling" << std::endl;
    file << "z_min " << M_z_min << std::endl;
    file << "z_max " << M_z_max << std::endl;
    file << "order " << M_polynomialOrder << std::endl;
    
    file << "pressure " << M_pressureCoefficients.size();
    for (const auto& c : M_pressureCoefficients) {
        file << " " << c;
    }
    file << std::endl;
    
    file << "disp_x " << M_displacementXCoefficients.size();
    for (const auto& c : M_displacementXCoefficients) {
        file << " " << c;
    }
    file << std::endl;
    
    file << "disp_y " << M_displacementYCoefficients.size();
    for (const auto& c : M_displacementYCoefficients) {
        file << " " << c;
    }
    file << std::endl;
    
    file << "disp_z " << M_displacementZCoefficients.size();
    for (const auto& c : M_displacementZCoefficients) {
        file << " " << c;
    }
    file << std::endl;
    
    file.close();
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

// ============================================================================
// Private Helper Methods
// ============================================================================
void RussianDollProblem::printCoefficients() const {
    std::cout << "[RussianDoll] Polynomial coefficients:" << std::endl;
    std::cout << "  Z-range: [" << M_z_min << ", " << M_z_max << "]" << std::endl;
    
    std::cout << "  Pressure: [";
    for (size_t i = 0; i < M_pressureCoefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << M_pressureCoefficients[i];
    }
    std::cout << "]" << std::endl;
    
    std::cout << "  Disp X: [";
    for (size_t i = 0; i < M_displacementXCoefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << M_displacementXCoefficients[i];
    }
    std::cout << "]" << std::endl;
    
    std::cout << "  Disp Y: [";
    for (size_t i = 0; i < M_displacementYCoefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << M_displacementYCoefficients[i];
    }
    std::cout << "]" << std::endl;
    
    std::cout << "  Disp Z: [";
    for (size_t i = 0; i < M_displacementZCoefficients.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << M_displacementZCoefficients[i];
    }
    std::cout << "]" << std::endl;
    
    // Print sample values
    std::cout << "  Sample p_v values:" << std::endl;
    std::cout << "    p(z_min=" << M_z_min << ") = " << evaluatePVPressure(M_z_min) << std::endl;
    std::cout << "    p(z_mid=" << (M_z_min+M_z_max)/2 << ") = " 
              << evaluatePVPressure((M_z_min+M_z_max)/2) << std::endl;
    std::cout << "    p(z_max=" << M_z_max << ") = " << evaluatePVPressure(M_z_max) << std::endl;
}