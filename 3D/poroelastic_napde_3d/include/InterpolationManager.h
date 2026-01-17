// ============================================================================
// InterpolationManager.h - Manages solution transfer between PV and PLC meshes
// ============================================================================
#ifndef INTERPOLATIONMANAGER_H
#define INTERPOLATIONMANAGER_H

#include "Core.h"
#include "Bulk.h"
#include <getfem/getfem_interpolation.h>
// Removed: #include <Eigen/Dense> - using GMM instead

/**
 * @class InterpolationManager
 * @brief Handles solution interpolation between PV and PLC meshes
 * 
 * Supports two approaches:
 * 1. LINE_INTERPOLATION: Extract along 1D lines, fit polynomial
 * 2. MESH_INTERPOLATION: Direct 3D-to-3D via GetFEM interpolation
 */

// Forward declarations
class Bulk;

// ============================================================================
// Line Profile Configuration
// ============================================================================
struct LineProfile {
    bgeot::base_node start_point;    ///< Line start (x,y,z)
    bgeot::base_node end_point;      ///< Line end (x,y,z)
    size_type num_samples;           ///< Number of sample points
    size_type polynomial_order;      ///< Order for OLS fit
    std::string name;                ///< Identifier for this line
    
    LineProfile() : num_samples(100), polynomial_order(3), name("line0") {}
};

// ============================================================================
// Polynomial Fit Result (GMM version - no Eigen)
// ============================================================================
struct PolynomialFit {
    std::vector<scalar_type> coefficients;    ///< Polynomial coefficients [a0, a1, ..., an]
    scalar_type r_squared;                     ///< R² goodness of fit
    scalar_type arc_length;                    ///< Total arc length of the line
    
    PolynomialFit() : r_squared(0.0), arc_length(0.0) {}
    
    /// Evaluate polynomial at normalized parameter t ∈ [0, 1]
    scalar_type evaluate(scalar_type t) const {
        scalar_type result = 0.0;
        scalar_type t_power = 1.0;
        for (size_type i = 0; i < coefficients.size(); ++i) {
            result += coefficients[i] * t_power;
            t_power *= t;
        }
        return result;
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
    
    // ========================================================================
    // Approach 1: LINE INTERPOLATION
    // ========================================================================
    
    /**
     * @brief Extract solution values along a line
     * @param solution FEM solution vector on source mesh
     * @param mf_source MeshFem for the solution
     * @param profile Line definition
     * @param[out] arc_coords Arc-length coordinates of samples
     * @param[out] values Solution values at sample points
     */
    void extractAlongLine(const scalarVectorPtr_Type& solution,
                          const getfem::mesh_fem& mf_source,
                          const LineProfile& profile,
                          std::vector<scalar_type>& arc_coords,
                          std::vector<scalar_type>& values);
    
    /**
     * @brief Fit polynomial to extracted data using OLS (GMM version)
     * @param arc_coords Arc-length coordinates (normalized to [0,1])
     * @param values Solution values
     * @param order Polynomial order
     * @return PolynomialFit result
     */
    PolynomialFit fitPolynomial(const std::vector<scalar_type>& arc_coords,
                                const std::vector<scalar_type>& values,
                                size_type order);
    
    /**
     * @brief Apply polynomial as boundary condition on target mesh
     * @param fit Polynomial fit result
     * @param profile Original line definition
     * @param mf_target Target mesh_fem
     * @param[out] bc_values Values at boundary DOFs
     */
    void applyPolynomialBC(const PolynomialFit& fit,
                           const LineProfile& profile,
                           const getfem::mesh_fem& mf_target,
                           scalarVectorPtr_Type& bc_values);
    
    /**
     * @brief Full line interpolation pipeline: extract → fit → apply
     * @param source_solution Solution on PV mesh
     * @param mf_source Source mesh_fem
     * @param target_solution Output: interpolated on PLC mesh/boundary
     * @param mf_target Target mesh_fem
     */
    void interpolateViaLines(const scalarVectorPtr_Type& source_solution,
                             const getfem::mesh_fem& mf_source,
                             scalarVectorPtr_Type& target_solution,
                             const getfem::mesh_fem& mf_target);
    
    // ========================================================================
    // Approach 2: MESH INTERPOLATION (3D to 3D)
    // ========================================================================
    
    /**
     * @brief Direct mesh-to-mesh interpolation using GetFEM
     * @param source_solution Solution on source (PV) mesh
     * @param mf_source Source mesh_fem
     * @param target_solution Output: interpolated on target (PLC) mesh
     * @param mf_target Target mesh_fem
     * 
     * Uses getfem::interpolation() internally
     */
    void interpolateViaMesh(const scalarVectorPtr_Type& source_solution,
                            const getfem::mesh_fem& mf_source,
                            scalarVectorPtr_Type& target_solution,
                            const getfem::mesh_fem& mf_target);
    
    /**
     * @brief Pre-compute interpolation matrix for efficiency
     * 
     * For repeated interpolations, compute once:
     * target_solution = M_interp_matrix * source_solution
     */
    void buildInterpolationMatrix(const getfem::mesh_fem& mf_source,
                                  const getfem::mesh_fem& mf_target);
    
    /**
     * @brief Fast interpolation using pre-computed matrix
     */
    void interpolateViaMatrix(const scalarVectorPtr_Type& source_solution,
                              scalarVectorPtr_Type& target_solution);
    
    // ========================================================================
    // Unified Interface
    // ========================================================================
    
    /**
     * @brief Interpolate from source to target using configured approach
     * @param source_solution Input solution
     * @param mf_source Source mesh_fem
     * @param target_solution Output solution
     * @param mf_target Target mesh_fem
     */
    void interpolate(const scalarVectorPtr_Type& source_solution,
                     const getfem::mesh_fem& mf_source,
                     scalarVectorPtr_Type& target_solution,
                     const getfem::mesh_fem& mf_target);
    
    // ========================================================================
    // Getters
    // ========================================================================
    
    inline Approach getApproach() const { return M_approach; }
    inline const std::vector<LineProfile>& getLineProfiles() const { return M_lineProfiles; }
    inline const std::vector<PolynomialFit>& getPolynomialFits() const { return M_polynomialFits; }
    
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
};

#endif // INTERPOLATIONMANAGER_H