// ============================================================================
// BC.h - Boundary condition management with polynomial callback support
// ============================================================================
#ifndef BC_H
#define BC_H

#include "Core.h"
#include "Parser.h"
#include <functional>

/**
 * @class BC
 * @brief Boundary condition handler for 3D bulk domain
 * 
 * Manages Dirichlet, Neumann, and mixed boundary conditions for
 * scalar and vector problems. Supports automatic boundary detection from
 * Gmsh physical tags or geometric criteria.
 * 
 * Enhanced with polynomial callback support for Russian Doll coupling:
 * - Pressure Neumann BC on outer wall: p_l = p_v(z)
 * - Displacement Dirichlet BC on outer wall: u_l = u_v(z)
 */
class BC {
public:
    // ========================================================================
    // Type Definitions for Callbacks
    // ========================================================================
    
    /**
     * @typedef ScalarCallback
     * @brief Function type for scalar BC as function of z
     */
    using ScalarCallback = std::function<scalar_type(scalar_type z)>;
    
    /**
     * @typedef VectorCallback  
     * @brief Function type for vector BC as function of z
     */
    using VectorCallback = std::function<bgeot::base_node(scalar_type z)>;

    // ========================================================================
    // Construction
    // ========================================================================
    
    /**
     * @brief Construct boundary conditions from data file
     * @param dataFile GetPot parameter file
     * @param problem Problem type identifier ("darcy/" or "mecc/")
     * @param section Data section for bulk domain (default: "bulkData/")
     */
    BC(const GetPot& dataFile,
       const std::string& problem,
       const std::string& section = "bulkData/");
    
    // ========================================================================
    // Standard BC Evaluation Functions
    // ========================================================================
    
    /// Evaluate scalar Neumann BC (natural boundary condition)
    scalar_type BCNeum(const base_node& x, const size_type& flag, const scalar_type t);
    
    /// Evaluate scalar Dirichlet BC (essential boundary condition)
    scalar_type BCDiri(const base_node& x, const size_type& flag);
    
    /// Evaluate vector-valued Dirichlet BC (e.g., for displacement)
    bgeot::base_node BCDiriVec(const base_node& x, const size_type& flag, const scalar_type t);
    
    /// Evaluate vector-valued velocity Dirichlet BC
    bgeot::base_node BCDiriVel(const base_node& x, const size_type& flag, const scalar_type t);
    
    /// Evaluate vector-valued Neumann BC (e.g., traction)
    bgeot::base_node BCNeumVec(const base_node& x, const size_type& flag, const scalar_type t);
    
    // ========================================================================
    // Polynomial Callback Methods (for Russian Doll Coupling)
    // ========================================================================
    
    /**
     * @brief Set pressure callback for a specific boundary region
     * @param region Boundary region ID (e.g., 0 for outer wall)
     * @param callback Function p(z) that returns pressure at height z
     * 
     * When set, BCNeum will use this callback instead of parser for the region
     */
    void setPressureCallback(size_type region, ScalarCallback callback);
    
    /**
     * @brief Set displacement callback for a specific boundary region
     * @param region Boundary region ID (e.g., 0 for outer wall)
     * @param callback Function u(z) that returns displacement vector at height z
     * 
     * When set, BCDiriVec will use this callback instead of parser for the region
     */
    void setDisplacementCallback(size_type region, VectorCallback callback);
    
    /**
     * @brief Set displacement callbacks from individual component functions
     * @param region Boundary region ID
     * @param callback_x Function for x-component u_x(z)
     * @param callback_y Function for y-component u_y(z)
     * @param callback_z Function for z-component u_z(z)
     */
    void setDisplacementCallbacks(size_type region,
                                  ScalarCallback callback_x,
                                  ScalarCallback callback_y,
                                  ScalarCallback callback_z);
    
    /**
     * @brief Set polynomial coefficients for pressure BC
     * @param region Boundary region ID
     * @param coefficients Polynomial coefficients [c0, c1, c2, ...]
     * @param z_min Minimum z for normalization
     * @param z_max Maximum z for normalization
     * 
     * Creates internal callback: p(z) = c0 + c1*t + c2*t^2 + ...
     * where t = (z - z_min) / (z_max - z_min)
     */
    void setPressurePolynomial(size_type region,
                               const std::vector<scalar_type>& coefficients,
                               scalar_type z_min,
                               scalar_type z_max);
    
    /**
     * @brief Set polynomial coefficients for displacement BC
     * @param region Boundary region ID
     * @param coeffs_x X-component coefficients
     * @param coeffs_y Y-component coefficients
     * @param coeffs_z Z-component coefficients
     * @param z_min Minimum z for normalization
     * @param z_max Maximum z for normalization
     */
    void setDisplacementPolynomial(size_type region,
                                   const std::vector<scalar_type>& coeffs_x,
                                   const std::vector<scalar_type>& coeffs_y,
                                   const std::vector<scalar_type>& coeffs_z,
                                   scalar_type z_min,
                                   scalar_type z_max);
    
    /**
     * @brief Clear pressure callback for a region
     * @param region Boundary region ID
     */
    void clearPressureCallback(size_type region);
    
    /**
     * @brief Clear displacement callback for a region
     * @param region Boundary region ID
     */
    void clearDisplacementCallback(size_type region);
    
    /**
     * @brief Check if pressure callback is set for a region
     */
    bool hasPressureCallback(size_type region) const;
    
    /**
     * @brief Check if displacement callback is set for a region
     */
    bool hasDisplacementCallback(size_type region) const;
    
    // ========================================================================
    // Boundary Region Management
    // ========================================================================
    
    /**
     * @brief Associate boundary regions with mesh (geometric detection)
     * @param meshPtr Pointer to mesh
     */
    void setBoundariesCylinder(getfem::mesh* meshPtr);
    
    void setBoundariesSquare(getfem::mesh* meshPtr);

    /**
     * @brief Associate boundary regions with mesh using Gmsh tags
     * @param meshPtr Pointer to mesh
     * @param regmap Map from physical names to region IDs (from Gmsh import)
     */
    void setBoundariesFromTagsName(getfem::mesh* meshPtr, 
                               const std::map<std::string, size_type>& regmap);
    
     /**
     * @brief Associate boundary regions with mesh using tag numbers in sorted order
     * 
     * This method selects the N largest (or smallest) tag numbers from the mesh
     * and assigns them to internal region IDs 0, 1, 2, ... in sorted order.
     * 
     * Example with 7 tags (1,2,3,4,5,6,7), 4 BCs, LARGEST mode:
     *   - Selects tags 4, 5, 6, 7
     *   - Internal region 0 <- Gmsh tag 4
     *   - Internal region 1 <- Gmsh tag 5
     *   - Internal region 2 <- Gmsh tag 6
     *   - Internal region 3 <- Gmsh tag 7
     * 
     * @param meshPtr Pointer to mesh
     * @param regmap Map from physical names to region IDs (from Gmsh import)
     * @param mode Whether to select LARGEST or SMALLEST tag numbers
     */
    void setBoundariesFromTagNumbers(getfem::mesh* meshPtr, 
                                     const std::map<std::string, size_type>& regmap,
                                     bool largest = true);
    
    /**
     * @brief Associate boundary regions with mesh using tag numbers (direct version)
     * 
     * This version works directly with the mesh regions without requiring the regmap.
     * It finds all boundary regions in the mesh and selects the N largest (or smallest).
     * 
     * @param meshPtr Pointer to mesh
     * @param mode Whether to select LARGEST or SMALLEST tag numbers
     */
    void setBoundariesFromTagNumbersDirect(getfem::mesh* meshPtr,
                                           bool largest=true );
    
    /**
     * @brief Get the mapping from internal region IDs to original Gmsh tag numbers
     * @return Map where key = internal region ID, value = original Gmsh tag
     * 
     * Useful for debugging and understanding which Gmsh tags were assigned.
     */
    const std::map<size_type, size_type>& getInternalToGmshMapping() const {
        return M_internalToGmsh;
    }
    
                               
    /// Get list of boundary region IDs with Neumann conditions
    std::vector<size_type> getNeumBD();
    
    /// Get list of boundary region IDs with Dirichlet conditions
    std::vector<size_type> getDiriBD();
    
    /// Get list of boundary region IDs with mixed conditions
    std::vector<size_type> getMixedBD();
    
    /// Get BC flag for a specific region (0=Dirichlet, 1=Neumann, 2=Mixed)
    size_type getBCFlag(size_type region) const {
        if (region < M_BC.size()) return M_BC[region];
        return 1; // Default to Neumann
    }

private:
    // ========================================================================
    // Helper Methods
    // ========================================================================
    
    /**
     * @brief Evaluate polynomial at normalized coordinate
     * @param z Physical z-coordinate
     * @param coeffs Polynomial coefficients
     * @param z_min Minimum z for normalization
     * @param z_max Maximum z for normalization
     * @return Polynomial value
     */
    scalar_type evaluatePolynomial(scalar_type z,
                                   const std::vector<scalar_type>& coeffs,
                                   scalar_type z_min,
                                   scalar_type z_max) const;

    /**
     * @brief Extract and sort tag numbers from regmap
     * @param regmap Map from physical names to region IDs
     * @return Sorted vector of unique tag numbers
     */
    std::vector<size_type> extractSortedTagNumbers(
        const std::map<std::string, size_type>& regmap) const;
    
    /**
     * @brief Select N tag numbers based on mode
     * @param sortedTags Sorted vector of all tag numbers
     * @param n Number of tags to select
     * @param mode Selection mode (LARGEST or SMALLEST)
     * @return Vector of selected tag numbers (sorted in ascending order)
     */
    std::vector<size_type> selectTagNumbers(
        const std::vector<size_type>& sortedTags,
        size_type n,
        bool largest) const;

    // ========================================================================
    // Member Variables
    // ========================================================================
    
    std::string M_section;      ///< Bulk data section name
    size_type M_nBoundaries;    ///< Number of boundary segments
    
    // Boundary condition expression strings
    std::string M_BCstring;     ///< BC descriptor (flags: 0=Dirichlet, 1=Neumann, 2=Mixed)
    std::string M_BCNeum;       ///< Scalar Neumann expression
    std::string M_BCNeumVec;    ///< Vector Neumann expression
    std::string M_BCDiri;       ///< Scalar Dirichlet expression
    std::string M_BCDiriVec;    ///< Vector Dirichlet expression
    std::string M_BCDiriVel;    ///< Velocity Dirichlet expression
    
    // Boundary region identifiers
    std::vector<size_type> M_NeumRG;    ///< Neumann region IDs
    std::vector<size_type> M_DiriRG;    ///< Dirichlet region IDs
    std::vector<size_type> M_MixedRG;   ///< Mixed region IDs
    std::vector<size_type> M_BC;        ///< All BC flags
    
    LifeV::Parser M_parser;     ///< Expression parser for BC evaluation
    
    // ========================================================================
    // Callback Storage (for polynomial BCs)
    // ========================================================================
    
    /// Pressure callbacks by region ID
    std::map<size_type, ScalarCallback> M_pressureCallbacks;
    
    /// Displacement callbacks by region ID (returns vector)
    std::map<size_type, VectorCallback> M_displacementCallbacks;
    
    /// Displacement component callbacks by region ID
    std::map<size_type, ScalarCallback> M_dispXCallbacks;
    std::map<size_type, ScalarCallback> M_dispYCallbacks;
    std::map<size_type, ScalarCallback> M_dispZCallbacks;
    
    /// Polynomial data storage (for regions using polynomial BC)
    struct PolynomialData {
        std::vector<scalar_type> coefficients;
        scalar_type z_min;
        scalar_type z_max;
    };
    
    std::map<size_type, PolynomialData> M_pressurePolynomials;
    std::map<size_type, PolynomialData> M_dispXPolynomials;
    std::map<size_type, PolynomialData> M_dispYPolynomials;
    std::map<size_type, PolynomialData> M_dispZPolynomials;

    
    /// Map from internal region ID to original Gmsh tag number
    std::map<size_type, size_type> M_internalToGmsh;
    
    /// Map from original Gmsh tag number to internal region ID
    std::map<size_type, size_type> M_gmshToInternal;
};

#endif // BC_H