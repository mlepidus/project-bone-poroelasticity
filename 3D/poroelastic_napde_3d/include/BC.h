// ============================================================================
// BC.h - Boundary condition management
// ============================================================================
#ifndef BC_H
#define BC_H

#include "Core.h"
#include "Parser.h"

/**
 * @class BC
 * @brief Boundary condition handler for bulk and fracture domains
 * 
 * Manages Dirichlet, Neumann, and mixed boundary conditions for both
 * scalar and vector problems. Generic implementation supporting multiple
 * physical problems.
 */
class BC {
public:
    /**
     * @brief Construct boundary conditions from data file
     * @param dataFile GetPot parameter file
     * @param problem Problem type identifier
     * @param section1 Data section for bulk domain
     * @param section2 Data section for fracture domain
     */
    BC(const GetPot& dataFile,
       const std::string& problem,
       const std::string& section1 = "bulkData/",
       const std::string& section2 = "fractureData/");
    
    /// Evaluate scalar Neumann BC (natural boundary condition)
    scalar_type BCNeum(const base_node& x, const size_type& flag, const scalar_type t);
    
    /// Evaluate scalar Neumann BC on fracture
    scalar_type BCNeumF(const base_node& x, const size_type& flag);
    
    /// Evaluate scalar Dirichlet BC (essential boundary condition)
    scalar_type BCDiri(const base_node& x, const size_type& flag);
    
    /// Evaluate scalar Dirichlet BC on fracture
    scalar_type BCDiriF(const base_node& x, const size_type& flag);
    
    /// Evaluate mixed BC (allows tangential slip, normal constraint)
    scalar_type BCMixed(const base_node& x, const size_type& flag);
    
    /// Evaluate vector-valued Dirichlet BC (e.g., for displacement)
    bgeot::base_node BCDiriVec(const base_node& x, const size_type& flag, const scalar_type t);
    
    /// Evaluate vector-valued velocity Dirichlet BC
    bgeot::base_node BCDiriVel(const base_node& x, const size_type& flag, const scalar_type t);
    
    /// Evaluate vector-valued Neumann BC (e.g., traction)
    bgeot::base_node BCNeumVec(const base_node& x, const size_type& flag, const scalar_type t);
    
    /**
     * @brief Associate boundary regions with mesh
     * @param meshPtr Pointer to mesh
     * @param where Domain selector: "bulk" or "fracture"
     */
    void setBoundaries(getfem::mesh* meshPtr, std::string where = "bulk");
    
    /**
     * @brief Associate boundary regions with mesh using Gmsh tags
     * @param meshPtr Pointer to mesh
     * @param regmap Map from physical names to region IDs (from Gmsh import)
     * @param where Domain selector: "bulk" or "fracture"
     */
    void setBoundariesFromTags(getfem::mesh* meshPtr, 
                               const std::map<std::string, size_type>& regmap,
                               std::string where = "bulk");
    


    /// Get list of boundary region IDs with Neumann conditions
    std::vector<size_type> getNeumBD(std::string where = "bulk");
    
    /// Get list of boundary region IDs with Dirichlet conditions
    std::vector<size_type> getDiriBD(std::string where = "bulk");
    
    /// Get list of boundary region IDs with mixed conditions
    std::vector<size_type> getMixedBD(std::string where = "bulk");

private:
    std::string M_section1;  ///< Bulk data section name
    std::string M_section2;  ///< Fracture data section name
    
    // Boundary condition expression strings
    std::string M_BCstring;      ///< Bulk BC descriptor
    std::string M_BCstringF;     ///< Fracture BC descriptor
    std::string M_BCNeum;        ///< Scalar Neumann expression
    std::string M_BCNeumF;       ///< Scalar Neumann expression (fracture)
    std::string M_BCNeumVec;     ///< Vector Neumann expression
    std::string M_BCDiri;        ///< Scalar Dirichlet expression
    std::string M_BCMixed;       ///< Mixed BC expression
    std::string M_BCDiriVec;     ///< Vector Dirichlet expression
    std::string M_BCDiriVel;     ///< Velocity Dirichlet expression
    std::string M_BCDiriF;       ///< Scalar Dirichlet expression (fracture)
    
    // Boundary region identifiers
    std::vector<size_type> M_NeumRG;   ///< Neumann region IDs (bulk)
    std::vector<size_type> M_DiriRG;   ///< Dirichlet region IDs (bulk)
    std::vector<size_type> M_MixedRG;  ///< Mixed region IDs (bulk)
    std::vector<size_type> M_NeumRGF;  ///< Neumann region IDs (fracture)
    std::vector<size_type> M_DiriRGF;  ///< Dirichlet region IDs (fracture)
    std::vector<size_type> M_BC;       ///< All BC region IDs (bulk)
    std::vector<size_type> M_BCF;      ///< All BC region IDs (fracture)
    
    LifeV::Parser M_parser;  ///< Expression parser for BC evaluation
    
    int M_nBoundaries;                   ///< Number of boundary segments
    std::vector<scalar_type> M_bdNodesX; ///< X-coordinates of boundary nodes
    std::vector<scalar_type> M_bdNodesY; ///< Y-coordinates of boundary nodes
};

#endif // BC_H