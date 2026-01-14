// ============================================================================
// Bulk.h - Bulk domain management
// ============================================================================
#ifndef BULK_H
#define BULK_H

#include "Core.h"
#include "Parser.h"
#include "BulkDarcyData.h"
#include "BulkElastData.h"
#include "BC.h"

/**
 * @class Bulk
 * @brief Manages 2D bulk domain including mesh and physical properties
 * 
 * Central container for mesh geometry, Darcy flow data, and elastic
 * material properties. Handles mesh construction and boundary regions.
 */
class Bulk {
public:
    /**
     * @brief Construct bulk domain from data file
     * @param dataFile GetPot parameter file
     * @param section Main bulk data section
     * @param sectionDomain Domain geometry section
     * @param sectionDarcy Darcy flow parameters section
     * @param sectionElast Elasticity parameters section
     */
    Bulk(const GetPot& dataFile,
         const std::string& section = "bulkData/",
         const std::string& sectionDomain = "domain/",
         const std::string& sectionDarcy = "darcy/",
         const std::string& sectionElast = "mecc/");
    
    /// Export mesh to VTK format for visualization
    void exportMesh(std::string filename);
    
    /// Get pointer to mesh object
    inline getfem::mesh* getMesh() { return &M_mesh; }
    
    /// Get vector of Neumann boundary regions
    inline std::vector<getfem::mesh_region*>* getNeumBoundaries() {
        return &M_NeumBoundaries;
    }
    
    /// Get pointer to Darcy problem data
    inline BulkDarcyData* getDarcyData() { return M_DarcyDataPtr; }
    
    /// Get pointer to elasticity problem data
    inline BulkElastData* getElastData() { return M_ElastDataPtr; }
    
    /// Get domain length in x-direction
    inline scalar_type Lx() { return M_Lx; }
    
    /// Get domain length in y-direction
    inline scalar_type Ly() { return M_Ly; }
    
    /// Get domain length in z-direction
    inline scalar_type Lz() { return M_Lz; }

    inline size_type getDim() {return M_dim;}
    /**
     * @brief Get the region map from Gmsh import
     * @return Map from physical names to region IDs
     */
    inline const std::map<std::string, size_type>& getRegionMap() const {
        return M_regmap;
    }
    
    /**
     * @brief Check if mesh was imported from external file
     * @return True if external mesh, false if generated internally
     */
    inline bool hasExternalMesh() const {
        return M_hasExternalMesh;
    }

private:
    // Data file sections
    std::string M_section;
    std::string M_sectionDomain;
    std::string M_sectionDarcy;
    std::string M_sectionElast;
    
    std::string M_meshFile;    ///< Mesh file name
    std::string M_meshFolder;  ///< Mesh file directory
    GetPot M_datafile;         ///< Parameter file handle
    
    BulkDarcyData M_DarcyData;      ///< Darcy flow data object
    BulkDarcyData* M_DarcyDataPtr;  ///< Pointer to Darcy data
    
    BulkElastData M_ElastData;      ///< Elasticity data object
    BulkElastData* M_ElastDataPtr;  ///< Pointer to elasticity data
    
    // Mesh discretization parameters
    size_type M_dim;              ///< dimension of the problem
    size_type M_Nx, M_Ny, M_Nz;  ///< Number of elements in each direction
    scalar_type M_Lx, M_Ly, M_Lz; ///< Domain dimensions
    
    // Nitsche's method penalty parameters
    scalar_type M_coeffNitscheNormal;   ///< Normal direction penalty
    scalar_type M_coeffNitscheTangent;  ///< Tangential direction penalty
    
    std::string M_meshType;   ///< Mesh type identifier
    getfem::mesh M_mesh;      ///< GetFem++ mesh object
    
    std::vector<getfem::mesh_region*> M_NeumBoundaries;  ///< Neumann boundary regions
    
    bool M_isVertical;  ///< Flag for vertical fracture orientation

    std::map<std::string, size_type> M_regmap;  ///< Physical names to region IDs
    bool M_hasExternalMesh;  ///< Flag indicating external mesh import
};

#endif // BULK_H
