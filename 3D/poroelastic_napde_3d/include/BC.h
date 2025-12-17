// ============================================================================
// BC.h - Boundary condition management
// ============================================================================
#ifndef BC_H
#define BC_H

#include "Core.h"
#include "Parser.h"

/**
 * @class BC
 * @brief Boundary condition handler for 3D bulk domain
 * 
 * Manages Dirichlet, Neumann, and mixed boundary conditions for
 * scalar and vector problems. Supports automatic boundary detection from
 * Gmsh physical tags or geometric criteria.
 */
class BC {
public:
    /**
     * @brief Construct boundary conditions from data file
     * @param dataFile GetPot parameter file
     * @param problem Problem type identifier ("darcy/" or "mecc/")
     * @param section Data section for bulk domain (default: "bulkData/")
     */
    BC(const GetPot& dataFile,
       const std::string& problem,
       const std::string& section = "bulkData/");
    
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
    
    /**
     * @brief Associate boundary regions with mesh (geometric detection)
     * @param meshPtr Pointer to mesh
     */
    void setBoundaries(getfem::mesh* meshPtr);
    
    /**
     * @brief Associate boundary regions with mesh using Gmsh tags
     * @param meshPtr Pointer to mesh
     * @param regmap Map from physical names to region IDs (from Gmsh import)
     */
    void setBoundariesFromTags(getfem::mesh* meshPtr, 
                               const std::map<std::string, size_type>& regmap);
    
    /// Get list of boundary region IDs with Neumann conditions
    std::vector<size_type> getNeumBD();
    
    /// Get list of boundary region IDs with Dirichlet conditions
    std::vector<size_type> getDiriBD();
    
    /// Get list of boundary region IDs with mixed conditions
    std::vector<size_type> getMixedBD();

private:
    std::string M_section;  ///< Bulk data section name
    
    // Boundary condition expression strings
    std::string M_BCstring;      ///< BC descriptor (flags: 0=Dirichlet, 1=Neumann, 2=Mixed)
    std::string M_BCNeum;        ///< Scalar Neumann expression
    std::string M_BCNeumVec;     ///< Vector Neumann expression
    std::string M_BCDiri;        ///< Scalar Dirichlet expression
    std::string M_BCDiriVec;     ///< Vector Dirichlet expression
    std::string M_BCDiriVel;     ///< Velocity Dirichlet expression
    
    // Boundary region identifiers
    std::vector<size_type> M_NeumRG;   ///< Neumann region IDs
    std::vector<size_type> M_DiriRG;   ///< Dirichlet region IDs
    std::vector<size_type> M_MixedRG;  ///< Mixed region IDs
    std::vector<size_type> M_BC;       ///< All BC flags
    
    LifeV::Parser M_parser;  ///< Expression parser for BC evaluation
    
    int M_nBoundaries;  ///< Number of boundary segments
};

#endif // BC_H