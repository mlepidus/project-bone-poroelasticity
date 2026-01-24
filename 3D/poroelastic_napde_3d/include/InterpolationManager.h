// ============================================================================
// InterpolationManager.h - Manages solution transfer between PV and PLC meshes
// ============================================================================
// Enhanced version with:
// - Line extraction for pressure and displacement
// - Polynomial fitting with explicit coefficient storage
// - Support for multi-field interpolation (pressure + displacement)
// - RHS coupling term computation
// - Boundary condition evaluation
// ============================================================================
#ifndef INTERPOLATIONMANAGER_H
#define INTERPOLATIONMANAGER_H

#include "Core.h"
#include "Bulk.h"
#include <getfem/getfem_interpolation.h>
#include <functional>

/**
 * @class InterpolationManager
 * @brief Handles solution interpolation between PV and PLC meshes
 * 
 * Supports two approaches:
 * 1. LINE_INTERPOLATION: Extract along 1D lines, fit polynomial, store coefficients
 * 2. MESH_INTERPOLATION: Direct 3D-to-3D via GetFEM interpolation
 * 
 * For Russian Doll coupling:
 * - Extracts PV pressure and displacement along a vertical line (z-axis)
 * - Fits polynomials: p_v(z) and u_v(z)
 * - Stores coefficients for BC evaluation on PLC outer wall
 * - Computes coupling terms for RHS modification
 */

// Forward declarations
class Bulk;

// ============================================================================
// Line Profile Configuration
// ============================================================================
/**
 * @struct LineProfile
 * @brief Configuration for extracting solution along a 1D line
 */
struct LineProfile {
    bgeot::base_node start_point;    ///< Line start (x,y,z)
    bgeot::base_node end_point;      ///< Line end (x,y,z)
    size_type num_samples;           ///< Number of sample points along line
    size_type polynomial_order;      ///< Order of polynomial to fit
    std::string name;                ///< Identifier for this line
    
    LineProfile() : num_samples(100), polynomial_order(3), name("line0") {
        start_point.resize(3);
        end_point.resize(3);
        for (int i = 0; i < 3; ++i) {
            start_point[i] = 0.0;
            end_point[i] = 0.0;
        }
        end_point[2] = 1.0;  // Default: vertical line from (0,0,0) to (0,0,1)
    }
};

// ============================================================================
// Polynomial Fit Result
// ============================================================================
/**
 * @struct PolynomialFit
 * @brief Result of polynomial fitting to extracted line data
 * 
 * The polynomial is expressed as: f(t) = c0 + c1*t + c2*t^2 + ... + cn*t^n
 * where t is the normalized parameter in [0, 1] along the line.
 */
struct PolynomialFit {
    std::vector<scalar_type> coefficients;    ///< Polynomial coefficients [c0, c1, ..., cn]
    scalar_type r_squared;                     ///< R² goodness of fit
    scalar_type arc_length;                    ///< Total arc length of the line
    scalar_type t_min;                         ///< Minimum parameter value (usually 0)
    scalar_type t_max;                         ///< Maximum parameter value (usually 1)
    scalar_type z_min;                         ///< Minimum z-coordinate (physical)
    scalar_type z_max;                         ///< Maximum z-coordinate (physical)
    std::string field_name;                    ///< Name of the field (pressure, disp_x, etc.)
    
    PolynomialFit() : r_squared(0.0), arc_length(0.0), 
                      t_min(0.0), t_max(1.0),
                      z_min(0.0), z_max(1.0),
                      field_name("unknown") {}
    
    /**
     * @brief Evaluate polynomial at normalized parameter t ∈ [0, 1]
     * @param t Normalized parameter
     * @return Polynomial value: c0 + c1*t + c2*t^2 + ...
     */
    scalar_type evaluate(scalar_type t) const {
        scalar_type result = 0.0;
        scalar_type t_power = 1.0;
        for (size_type i = 0; i < coefficients.size(); ++i) {
            result += coefficients[i] * t_power;
            t_power *= t;
        }
        return result;
    }
    
    /**
     * @brief Evaluate polynomial at physical coordinate z
     * @param z Physical z-coordinate
     * @return Polynomial value at normalized position
     */
    scalar_type evaluateAtZ(scalar_type z) const {
        scalar_type t = 0.0;
        if (std::abs(z_max - z_min) > 1e-15) {
            t = (z - z_min) / (z_max - z_min);
        }
        t = std::max(0.0, std::min(1.0, t));  // Clamp
        return evaluate(t);
    }
    
    /**
     * @brief Evaluate polynomial derivative at normalized parameter t
     * @param t Normalized parameter
     * @return Derivative value: c1 + 2*c2*t + 3*c3*t^2 + ...
     */
    scalar_type evaluateDerivative(scalar_type t) const {
        scalar_type result = 0.0;
        scalar_type t_power = 1.0;
        for (size_type i = 1; i < coefficients.size(); ++i) {
            result += i * coefficients[i] * t_power;
            t_power *= t;
        }
        return result;
    }
};

// ============================================================================
// Multi-Field Polynomial Storage
// ============================================================================
/**
 * @struct FieldPolynomials
 * @brief Stores polynomial fits for all fields (pressure + 3 displacement components)
 */
struct FieldPolynomials {
    PolynomialFit pressure;          ///< Pressure polynomial p_v(z)
    PolynomialFit displacement_x;    ///< Displacement x-component u_x(z)
    PolynomialFit displacement_y;    ///< Displacement y-component u_y(z)
    PolynomialFit displacement_z;    ///< Displacement z-component u_z(z)
    
    scalar_type z_min;               ///< Common z-range minimum
    scalar_type z_max;               ///< Common z-range maximum
    
    FieldPolynomials() : z_min(0.0), z_max(1.0) {}
    
    /**
     * @brief Evaluate all fields at a given z
     * @param z Physical z-coordinate
     * @param[out] p Pressure value
     * @param[out] ux Displacement x
     * @param[out] uy Displacement y
     * @param[out] uz Displacement z
     */
    void evaluateAll(scalar_type z, 
                     scalar_type& p, 
                     scalar_type& ux, 
                     scalar_type& uy, 
                     scalar_type& uz) const {
        p = pressure.evaluateAtZ(z);
        ux = displacement_x.evaluateAtZ(z);
        uy = displacement_y.evaluateAtZ(z);
        uz = displacement_z.evaluateAtZ(z);
    }
};

// ============================================================================
// InterpolationManager Class
// ============================================================================
class InterpolationManager {
public:
    enum class Approach {
        LINE_INTERPOLATION,
        MESH_INTERPOLATION
    };
    
    /**
     * @brief Constructor
     * @param dataFile GetPot parameter file
     * @param sourceBulk Source mesh domain (PV)
     * @param targetBulk Target mesh domain (PLC)
     */
    InterpolationManager(const GetPot& dataFile,
                         Bulk* sourceBulk,
                         Bulk* targetBulk);
    
    ~InterpolationManager() = default;
    
    // ========================================================================
    // Setup
    // ========================================================================
    
    /**
     * @brief Initialize interpolation structures
     * Reads line definitions from data file if using LINE approach
     */
    void initialize();
    
    /**
     * @brief Add a line profile for extraction (LINE approach)
     * @param profile Line configuration
     */
    void addLineProfile(const LineProfile& profile);
    
    /**
     * @brief Clear all line profiles
     */
    void clearLineProfiles() { M_lineProfiles.clear(); }
    
    // ========================================================================
    // Core Line Extraction Methods
    // ========================================================================
    
    /**
     * @brief Extract scalar solution values along a line
     * @param solution FEM solution vector on source mesh
     * @param mf_source MeshFem for the solution
     * @param profile Line definition
     * @param[out] arc_coords Arc-length coordinates of samples (normalized to [0,1])
     * @param[out] values Solution values at sample points
     */
    void extractAlongLine(const scalarVectorPtr_Type& solution,
                          const getfem::mesh_fem& mf_source,
                          const LineProfile& profile,
                          std::vector<scalar_type>& arc_coords,
                          std::vector<scalar_type>& values);
    
    /**
     * @brief Extract vector solution (displacement) along a line
     * @param solution FEM solution vector (interleaved ux,uy,uz)
     * @param mf_source MeshFem for the solution (qdim=3)
     * @param profile Line definition
     * @param[out] arc_coords Arc-length coordinates
     * @param[out] values_x X-component values
     * @param[out] values_y Y-component values
     * @param[out] values_z Z-component values
     */
    void extractVectorAlongLine(const scalarVectorPtr_Type& solution,
                                const getfem::mesh_fem& mf_source,
                                const LineProfile& profile,
                                std::vector<scalar_type>& arc_coords,
                                std::vector<scalar_type>& values_x,
                                std::vector<scalar_type>& values_y,
                                std::vector<scalar_type>& values_z);
    
    // ========================================================================
    // Polynomial Fitting Methods
    // ========================================================================
    
    /**
     * @brief Fit polynomial to extracted data using OLS (Ordinary Least Squares)
     * @param arc_coords Arc-length coordinates (normalized to [0,1])
     * @param values Solution values
     * @param order Polynomial order
     * @param field_name Name of the field being fitted
     * @return PolynomialFit result with coefficients
     */
    PolynomialFit fitPolynomial(const std::vector<scalar_type>& arc_coords,
                                const std::vector<scalar_type>& values,
                                size_type order,
                                const std::string& field_name = "unknown");
    
    /**
     * @brief Fit polynomial with specified z-range
     * @param z_coords Physical z-coordinates
     * @param values Solution values
     * @param order Polynomial order
     * @param z_min Minimum z for normalization
     * @param z_max Maximum z for normalization
     * @param field_name Name of the field
     * @return PolynomialFit result
     */
    PolynomialFit fitPolynomialWithZRange(const std::vector<scalar_type>& z_coords,
                                          const std::vector<scalar_type>& values,
                                          size_type order,
                                          scalar_type z_min,
                                          scalar_type z_max,
                                          const std::string& field_name = "unknown");
    
    // ========================================================================
    // Complete Interpolation Pipeline
    // ========================================================================
    
    /**
     * @brief Extract pressure and fit polynomial
     * @param pressure_solution Pressure solution from PV
     * @param mf_pressure Pressure mesh_fem
     * @return PolynomialFit for pressure
     */
    PolynomialFit interpolatePressure(const scalarVectorPtr_Type& pressure_solution,
                                      const getfem::mesh_fem& mf_pressure);
    
    /**
     * @brief Extract displacement and fit polynomials for all 3 components
     * @param displacement_solution Displacement solution from PV
     * @param mf_displacement Displacement mesh_fem
     * @param[out] fit_x Polynomial fit for x-component
     * @param[out] fit_y Polynomial fit for y-component
     * @param[out] fit_z Polynomial fit for z-component
     */
    void interpolateDisplacement(const scalarVectorPtr_Type& displacement_solution,
                                 const getfem::mesh_fem& mf_displacement,
                                 PolynomialFit& fit_x,
                                 PolynomialFit& fit_y,
                                 PolynomialFit& fit_z);
    
    /**
     * @brief Full interpolation of pressure and displacement
     * @param pressure_solution Pressure from PV
     * @param mf_pressure Pressure mesh_fem
     * @param displacement_solution Displacement from PV
     * @param mf_displacement Displacement mesh_fem
     * @return FieldPolynomials structure with all fits
     */
    FieldPolynomials interpolateAllFields(const scalarVectorPtr_Type& pressure_solution,
                                          const getfem::mesh_fem& mf_pressure,
                                          const scalarVectorPtr_Type& displacement_solution,
                                          const getfem::mesh_fem& mf_displacement);
    
    // ========================================================================
    // BC Application Methods
    // ========================================================================
    
    /**
     * @brief Apply polynomial as boundary condition on target mesh
     * @param fit Polynomial fit result
     * @param profile Original line definition
     * @param mf_target Target mesh_fem
     * @param[out] bc_values Values at all DOFs (based on z-projection)
     */
    void applyPolynomialBC(const PolynomialFit& fit,
                           const LineProfile& profile,
                           const getfem::mesh_fem& mf_target,
                           scalarVectorPtr_Type& bc_values);
    
    /**
     * @brief Create a BC callback function from polynomial fit
     * @param fit Polynomial fit
     * @return Lambda that evaluates p_v(z)
     */
    std::function<scalar_type(scalar_type)> createBCCallback(const PolynomialFit& fit);
    
    // ========================================================================
    // Legacy Methods for Mesh Interpolation
    // ========================================================================
    
    void interpolateViaLines(const scalarVectorPtr_Type& source_solution,
                             const getfem::mesh_fem& mf_source,
                             scalarVectorPtr_Type& target_solution,
                             const getfem::mesh_fem& mf_target);
    
    void interpolateViaMesh(const scalarVectorPtr_Type& source_solution,
                            const getfem::mesh_fem& mf_source,
                            scalarVectorPtr_Type& target_solution,
                            const getfem::mesh_fem& mf_target);
    
    void buildInterpolationMatrix(const getfem::mesh_fem& mf_source,
                                  const getfem::mesh_fem& mf_target);
    
    void interpolateViaMatrix(const scalarVectorPtr_Type& source_solution,
                              scalarVectorPtr_Type& target_solution);
    
    void interpolate(const scalarVectorPtr_Type& source_solution,
                     const getfem::mesh_fem& mf_source,
                     scalarVectorPtr_Type& target_solution,
                     const getfem::mesh_fem& mf_target);
    
    // ========================================================================
    // Coefficient Access (for BC evaluation in RussianDoll)
    // ========================================================================
    
    /**
     * @brief Get the last polynomial fit result (pressure)
     * @return Const reference to the most recent pressure polynomial fit
     */
    const PolynomialFit& getLastPressureFit() const;
    
    /**
     * @brief Get all field polynomials from last interpolation
     * @return Const reference to FieldPolynomials structure
     */
    const FieldPolynomials& getFieldPolynomials() const { return M_fieldPolynomials; }
    
    /**
     * @brief Get polynomial coefficients from the last pressure fit
     * @return Const reference to coefficient vector [c0, c1, ...]
     */
    const std::vector<scalar_type>& getPressureCoefficients() const;
    
    /**
     * @brief Get displacement coefficients (x-component)
     */
    const std::vector<scalar_type>& getDisplacementXCoefficients() const;
    
    /**
     * @brief Get displacement coefficients (y-component)
     */
    const std::vector<scalar_type>& getDisplacementYCoefficients() const;
    
    /**
     * @brief Get displacement coefficients (z-component)
     */
    const std::vector<scalar_type>& getDisplacementZCoefficients() const;
    
    /**
     * @brief Get z-range from the last line extraction
     * @param[out] z_min Minimum z value
     * @param[out] z_max Maximum z value
     */
    void getZRange(scalar_type& z_min, scalar_type& z_max) const;
    
    /**
     * @brief Evaluate the fitted pressure polynomial at a given z
     * @param z Physical z-coordinate
     * @return Interpolated pressure value
     */
    scalar_type evaluatePressureAtZ(scalar_type z) const;
    
    /**
     * @brief Evaluate the fitted displacement at a given z
     * @param z Physical z-coordinate
     * @param[out] ux X-component
     * @param[out] uy Y-component
     * @param[out] uz Z-component
     */
    void evaluateDisplacementAtZ(scalar_type z, 
                                 scalar_type& ux, 
                                 scalar_type& uy, 
                                 scalar_type& uz) const;
    
    // ========================================================================
    // Getters
    // ========================================================================
    
    inline Approach getApproach() const { return M_approach; }
    inline const std::vector<LineProfile>& getLineProfiles() const { return M_lineProfiles; }
    inline const std::vector<PolynomialFit>& getPolynomialFits() const { return M_polynomialFits; }
    inline bool isMatrixBuilt() const { return M_matrixBuilt; }
    
    // Legacy compatibility
    const PolynomialFit& getLastPolynomialFit() const;
    const std::vector<scalar_type>& getPolynomialCoefficients() const;
    scalar_type evaluateAtZ(scalar_type z) const;
    
private:
    // Configuration
    Approach M_approach;
    std::string M_section;
    
    // Domain pointers
    Bulk* M_sourceBulk;      ///< PV (outer) mesh
    Bulk* M_targetBulk;      ///< PLC (inner) mesh
    
    // Line interpolation data
    std::vector<LineProfile> M_lineProfiles;
    std::vector<PolynomialFit> M_polynomialFits;
    
    // Multi-field storage
    FieldPolynomials M_fieldPolynomials;
    
    // Z-range from last extraction
    scalar_type M_z_min;
    scalar_type M_z_max;
    
    // Mesh interpolation data
    gmm::row_matrix<gmm::wsvector<scalar_type>> M_interpMatrix;
    bool M_matrixBuilt;
    
    // Helper methods
    void readLineProfilesFromFile(const GetPot& dataFile);
    bgeot::base_node interpolatePoint(const LineProfile& profile, scalar_type t);
    
    // GMM helper for solving normal equations (replaces Eigen QR)
    void solveNormalEquations(const gmm::dense_matrix<scalar_type>& VtV,
                              const std::vector<scalar_type>& Vty,
                              std::vector<scalar_type>& coeffs);
    
    // Compute z-coordinates from sample points
    void computeZCoordinates(const std::vector<bgeot::base_node>& sample_points,
                             std::vector<scalar_type>& z_coords);
};

#endif // INTERPOLATIONMANAGER_H