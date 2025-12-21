#include "../include/Bulk.h"

Bulk::Bulk ( const GetPot& dataFile,
                             const std::string& section,
                             const std::string& sectionDomain,
                             const std::string& sectionDarcy,
                             const std::string& sectionElast
                             ) :
            M_section ( section ),
            M_sectionDomain ( M_section + sectionDomain ),
            M_sectionDarcy ( M_section + sectionDarcy ),
            M_sectionElast ( M_section + sectionElast ),   
	        M_meshFile( dataFile ( ( M_sectionDomain + "meshExternal" ).data (), "none" )   ),            
	        M_meshFolder( dataFile ( ( M_sectionDomain + "meshFolder" ).data (), "" )   ),

            M_datafile(dataFile),    
    	    M_DarcyData( dataFile), 
            M_DarcyDataPtr(nullptr),
	        M_ElastData( dataFile),
            M_ElastDataPtr(nullptr),


        M_Nx ( dataFile ( ( M_sectionDomain + "spatialDiscretizationX" ).data (), 10 ) ),
        M_Ny ( dataFile ( ( M_sectionDomain + "spatialDiscretizationY" ).data (), 10 ) ),
        M_Nz ( dataFile ( ( M_sectionDomain + "spatialDiscretizationZ" ).data (), 10 ) ),
	    M_Lx ( dataFile ( ( M_sectionDomain + "lengthAbscissa" ).data (), 1. ) ),
        M_Ly ( dataFile ( ( M_sectionDomain + "lengthOrdinate" ).data (), 1. ) ),
	    M_Lz ( dataFile ( ( M_sectionDomain + "lengthQuota" ).data (), 1. ) ),
	    

	    M_coeffNitscheNormal ( dataFile ( ( M_sectionDomain + "coeffN" ).data (), 1. ) ),
	    M_coeffNitscheTangent ( dataFile ( ( M_sectionDomain + "coeffT" ).data (), 1. ) ),
        M_meshType( dataFile ( ( M_sectionDomain + "meshType" ).data (), "GT_PK(3,1)" ) ),

	    M_isVertical( dataFile ( ( M_sectionDomain + "isVertical" ).data (), true ) ),
  
    
        M_hasExternalMesh(false)
{  

    bgeot::pgeometric_trans pgt; 

    pgt =  bgeot::geometric_trans_descriptor(M_meshType);
    size_type N = pgt->dim();

    std::vector<size_type> nsubdiv(N);
    
    // Check if external mesh file is specified
    if (M_meshFile.compare("none")==0)
    {
    // Generate internal structured mesh
    std::cout << "Creating a semi-structured default mesh "<<std::endl;
    nsubdiv[0]= M_Nx; 
    nsubdiv[1]= M_Ny;
    nsubdiv[2]= M_Nz;
    getfem::regular_unit_mesh(M_mesh, nsubdiv, pgt,false);  //creates a semi-structured mesh

    bgeot::base_matrix M(N,N);  //transformation matrix (scaling, shearing, etc)
    M(0,0)=M_Lx;  	        // scale the unit mesh to [LX,LY]
    M(1,1)=M_Ly; 
    M(2,2)=M_Lz; 
  
    M_mesh.transformation(M);
    M_hasExternalMesh = false;

    }
    else
    {
    // Import mesh from Gmsh file
       std::cout << "Importing mesh from file "<<std::endl;
	   std::cout << M_meshFile<<std::endl;
       //typedef std::map<std::string, size_type> RegMap; //to be checked
       //typedef RegMap::iterator RegMapIter;
       //RegMap regmap;
       getfem::import_mesh_gmsh(M_meshFolder + M_meshFile, M_mesh, M_regmap); //read mesh and boundary regions
       M_hasExternalMesh = true;
       //std::cout << M_regmap.size() << "\n";
        std::cout << "Mesh import successful!" << std::endl;
        std::cout << "Number of physical groups: " << M_regmap.size() << std::endl;
        
        // ====================================================================
        // Display imported physical names
        // ====================================================================
        std::cout << "\n--- Physical Names from Gmsh ---" << std::endl;
        for (const auto& pair : M_regmap) {
            std::cout << "  '" << pair.first << "' -> region ID " << pair.second << std::endl;
        }
        
        // ====================================================================
        // Verify regions exist in mesh
        // ====================================================================
        std::cout << "\n--- Verifying Mesh Regions ---" << std::endl;
        dal::bit_vector mesh_regions = M_mesh.regions_index();
        
        if (mesh_regions.empty()) {
            std::cout << "  WARNING: No regions found in mesh!" << std::endl;
        } else {
            std::cout << "Active regions in mesh:" << std::endl;
            for (dal::bv_visitor i(mesh_regions); !i.finished(); ++i) {
                // Count elements in this region
                getfem::mesh_region region = M_mesh.region(i);
                size_t count = 0;
                for (getfem::mr_visitor it(region); !it.finished(); ++it) {
                    count++;
                }
                std::cout << "  Region " << i << ": " << count << " elements" << std::endl;
            }
        }
        
        // ====================================================================
        // Check consistency between regmap and actual regions
        // ====================================================================
        std::cout << "\n--- Checking regmap consistency ---" << std::endl;
        for (const auto& [name, id] : M_regmap) {
            if (mesh_regions.is_in(id)) {
                std::cout << "  ✓ '" << name << "' (region " << id << ") EXISTS" << std::endl;
            } else {
                std::cout << "  ✗ WARNING: '" << name << "' (region " << id 
                          << ") NOT FOUND in mesh!" << std::endl;
            }
        }
        
        std::cout << "================================================\n" << std::endl;
    }

    M_DarcyDataPtr = &M_DarcyData;
    M_ElastDataPtr = &M_ElastData;
}

void Bulk::exportMesh(std::string nomefile)
{
   getfem::vtk_export vtkmesh(nomefile);
   vtkmesh.exporting(M_mesh);
   vtkmesh.write_mesh();

}




