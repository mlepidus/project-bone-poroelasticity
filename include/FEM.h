// ============================================================================
// FEM.h - Finite element method wrapper class
// ============================================================================
#ifndef FEM_H
#define FEM_H

#include "Core.h"

/**
 * @class FEM
 * @brief Wrapper class for GetFem++ mesh_fem providing simplified interface
 * 
 * This class encapsulates GetFem++'s finite element space functionality,
 * providing easier access to DOF information, spatial coordinates, and
 * extended DOF tracking for XFEM/level set applications.
 */
class FEM {
public:
    /**
     * @brief Construct FEM from data file configuration
     * @param mesh Pointer to the underlying mesh
     * @param dataFile GetPot parameter file
     * @param problem Problem identifier string
     * @param variable Variable name (e.g., "Pressure", "Displacement")
     * @param section Data file section containing FEM parameters
     */
    FEM(const getfem::mesh* mesh,
        const GetPot& dataFile,
        const std::string& problem,
        const std::string& variable,
        const std::string& section = "bulkData/");
    
    /**
     * @brief Construct FEM with explicit type specification
     * @param mesh Pointer to the underlying mesh
     * @param femType Finite element type string (e.g., "FEM_PK(2,1)")
     * @param spaceDim Spatial dimension of the vector field
     */
    FEM(const getfem::mesh* mesh,
        std::string femType,
        size_type spaceDim);

    /**
     * @brief Get number of degrees of freedom
     * @param which Selector: "all" (base + extended), "base", or "extended"
     * @return Total number of DOFs
     */
    size_type nb_dof(std::string which = "all");

    /// Get finite element type string
    inline std::string type() { return M_femType; }
    
    /// Get pointer to underlying GetFem++ mesh_fem object
    inline getfem::mesh_fem* getFEM() { return &M_FEM; }
    
    /// Get spatial coordinates of all DOF points
    inline std::vector<base_node> getDOFpoints() { return M_DOFpoints; }
    
    /// Get indices of extended DOFs (for XFEM enrichment)
    inline std::vector<size_type> getExt() { return M_extended; }
    
    /// Get signs (+1/-1) indicating which side of level set each extended DOF belongs to
    inline std::vector<scalar_type> getExtSign() { return M_extSign; }
    
    /// Get base DOF index corresponding to i-th extended DOF
    inline size_type getExt(size_type i) { return M_extended[i]; }
    
    /// Get sign of i-th DOF relative to level set
    inline scalar_type getDOFSign(size_type i) { return M_DOFSign[i]; }
    
    /// Get sign of i-th extended DOF
    inline scalar_type getExtSign(size_type i) { return M_extSign[i]; }
    
    /// Get spatial coordinate of i-th basic DOF
    inline base_node point_of_basic_dof(size_type i) {
        return M_FEM.point_of_basic_dof(i);
    }

    /// Get target dimension of the finite element space (scalar=1, vector=dim)
    size_type get_qdim() const { return M_FEM.get_qdim(); }

private:
    std::string M_section;              ///< Data file section name
    std::string M_femType;              ///< Finite element type descriptor
    size_type M_SpaceDim;               ///< Spatial dimension
    getfem::mesh_fem M_FEM;             ///< Underlying GetFem++ FEM object
    const getfem::mesh* M_meshPtr;      ///< Pointer to mesh
    std::vector<base_node> M_DOFpoints; ///< Spatial coordinates of DOF points
    
    // XFEM/enrichment data structures
    std::vector<size_type> M_extended;   ///< Indices of extended DOFs
    std::vector<scalar_type> M_extSign;  ///< Signs of extended DOFs
    std::vector<scalar_type> M_DOFSign;  ///< Signs of all DOFs
};

#endif // FEM_H