#include "../include/Bulk.h"

Bulk::Bulk ( const GetPot& dataFile,
                             const std::string& section,
                             const std::string& sectionDomain,
                             const std::string& sectionDarcy,
                             const std::string& sectionElast
                             ) :
            M_datafile(dataFile),                 
            M_section ( section ),
            M_sectionDomain ( M_section + sectionDomain ),
            M_sectionDarcy ( M_section + sectionDarcy ),
            M_sectionElast ( M_section + sectionElast ),
	    M_DarcyData( dataFile),
	    M_ElastData( dataFile),
            // domain
	    M_meshType( dataFile ( ( M_sectionDomain + "meshType" ).data (), "GT_PK(3,1)" ) ),
	    M_meshFile( dataFile ( ( M_sectionDomain + "meshExternal" ).data (), "none" )   ),
	    M_meshFolder( dataFile ( ( M_sectionDomain + "meshFolder" ).data (), "" )   ),
        M_Nx ( dataFile ( ( M_sectionDomain + "spatialDiscretizationX" ).data (), 10 ) ),
        M_Ny ( dataFile ( ( M_sectionDomain + "spatialDiscretizationY" ).data (), 10 ) ),
        M_Nz ( dataFile ( ( M_sectionDomain + "spatialDiscretizationZ" ).data (), 10 ) ),
	    M_Lx ( dataFile ( ( M_sectionDomain + "lengthAbscissa" ).data (), 1. ) ),
        M_Ly ( dataFile ( ( M_sectionDomain + "lengthOrdinate" ).data (), 1. ) ),
	    M_Lz ( dataFile ( ( M_sectionDomain + "lengthQuota" ).data (), 1. ) ),
	    M_coeffNitscheNormal ( dataFile ( ( M_sectionDomain + "coeffN" ).data (), 1. ) ),
	    M_coeffNitscheTangent ( dataFile ( ( M_sectionDomain + "coeffT" ).data (), 1. ) ),
	    M_isVertical( dataFile ( ( M_sectionDomain + "isVertical" ).data (), true ) )
{  

    bgeot::pgeometric_trans pgt; 

    pgt =  bgeot::geometric_trans_descriptor(M_meshType);
    size_type N = pgt->dim();

    std::vector<size_type> nsubdiv(N);
    
    if (M_meshFile.compare("none")==0)
    {
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
    }
    else
    {
      std::cout << "Importing mesh from file "<<std::endl;
	   std::cout << M_meshFile<<std::endl;
      // getfem::import_mesh(M_meshFolder + M_meshFile, "gmsh", M_mesh);
       typedef std::map<std::string, size_type> RegMap;
       typedef RegMap::iterator RegMapIter;
       RegMap regmap;
       getfem::import_mesh_gmsh(M_meshFolder + M_meshFile, M_mesh, regmap); //read mesh and boundary regions
       std::cout << regmap.size() << "\n";

       
// 1. What's in regmap
std::cout << "In regmap (from PhysicalNames):" << std::endl;
for (const auto& pair : regmap) {
    std::cout << "  '" << pair.first << "' -> id " << pair.second << std::endl;
}

// 2. What regions actually exist in mesh
std::cout << "\nActual regions in mesh (regions_index()):" << std::endl;
dal::bit_vector mesh_regions = M_mesh.regions_index();
if (mesh_regions.empty()) {
    std::cout << "  NO REGIONS IN MESH!" << std::endl;
} else {
    for (dal::bv_visitor i(mesh_regions); !i.finished(); ++i) {
        std::cout << "  Region id: " << i;
        
        // Count elements in this region
        getfem::mesh_region region = M_mesh.region(i);
        size_t count = 0;
        for (getfem::mr_visitor it(region); !it.finished(); ++it) {
            count++;
        }
        std::cout << " (" << count << " elements)" << std::endl;
    }
}

// 3. Check if regmap IDs exist in mesh
std::cout << "\nChecking regmap IDs against mesh:" << std::endl;
for (const auto& [name, id] : regmap) {
    if (mesh_regions.is_in(id)) {
        std::cout << "  ✓ '" << name << "' (id=" << id << ") EXISTS in mesh" << std::endl;
    } else {
        std::cout << "  ✗ '" << name << "' (id=" << id << ") NOT FOUND in mesh!" << std::endl;
    }
}
}


    M_DarcyDataPtr=&M_DarcyData;
    M_ElastDataPtr=&M_ElastData;

    

    
}

void Bulk::exportMesh(std::string nomefile)
{
   getfem::vtk_export vtkmesh(nomefile);
   vtkmesh.exporting(M_mesh);
   vtkmesh.write_mesh();

}




