#ifndef BULK_H
#define BULK_H 

#include "Core.h"
#include "Parser.h"
#include "BulkDarcyData.h"
#include "BulkElastData.h"

#include "BC.h"

//gestione del dominio 2D "bulk"

class Bulk
{
public: 
	 Bulk ( const GetPot& dataFile,
                   const std::string& section = "bulkData/",
                   const std::string& sectionDomain = "domain/",
                   const std::string& sectionDarcy = "darcy/",
                   const std::string& sectionElast = "mecc/");

  void exportMesh(std::string nomefile);
  
  inline getfem::mesh* getMesh()
  {
	return &M_mesh;
  }


  
  inline std::vector<getfem::mesh_region*>* getNeumBoundaries()
  {
	return &M_NeumBoundaries;
  }
  
  inline BulkDarcyData* getDarcyData()
  {
	return M_DarcyDataPtr;
  }
  
  inline BulkElastData* getElastData()
  {
	return M_ElastDataPtr;
  } 
  
  inline scalar_type Lx()
  {
  	return M_Lx;
  }
  
   inline scalar_type Ly()
  {
  	return M_Ly;
  }



private:

    // Attributes
    std::string M_section;
    std::string M_sectionDomain;
    std::string M_sectionDarcy;
    std::string M_sectionElast;
    std::string M_meshFile;
    std::string M_meshFolder;
    GetPot M_datafile;
    BulkDarcyData M_DarcyData;
    BulkDarcyData* M_DarcyDataPtr;
    
    BulkElastData M_ElastData;
    BulkElastData* M_ElastDataPtr;

    size_type M_Nx;
    size_type M_Ny;
    scalar_type M_Lx;
    scalar_type M_Ly;

    scalar_type M_coeffNitscheNormal;
    scalar_type M_coeffNitscheTangent;

    std::string M_meshType;
    getfem::mesh M_mesh;          // the mesh
    
    std::vector<getfem::mesh_region*> M_NeumBoundaries;
    
    bool M_isVertical;


};

#endif
