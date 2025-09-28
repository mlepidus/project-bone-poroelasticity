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
	    M_meshType( dataFile ( ( M_sectionDomain + "meshType" ).data (), "GT_PK(2,1)" ) ),
	    M_meshFile( dataFile ( ( M_sectionDomain + "meshExternal" ).data (), "none" )   ),
	    M_meshFolder( dataFile ( ( M_sectionDomain + "meshFolder" ).data (), "" )   ),
            M_Nx ( dataFile ( ( M_sectionDomain + "spatialDiscretizationX" ).data (), 10 ) ),
            M_Ny ( dataFile ( ( M_sectionDomain + "spatialDiscretizationY" ).data (), 10 ) ),
	    M_Lx ( dataFile ( ( M_sectionDomain + "lengthAbscissa" ).data (), 1. ) ),
            M_Ly ( dataFile ( ( M_sectionDomain + "lengthOrdinate" ).data (), 1. ) ),
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
    nsubdiv[0]= M_Nx; 
    nsubdiv[1]= M_Ny;
    getfem::regular_unit_mesh(M_mesh, nsubdiv, pgt,false);  //creates a semi-structured mesh

    bgeot::base_matrix M(N,N);  //transformation matrix (scaling, shearing, etc)
    M(0,0)=M_Lx;  	        // scale the unit mesh to [LX,LY]
    M(1,1)=M_Ly; 
  
    M_mesh.transformation(M);
    }
    else
    {
	std::cout << M_meshFile<<std::endl;
       getfem::import_mesh(M_meshFolder + M_meshFile,"gmsh", M_mesh);
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




