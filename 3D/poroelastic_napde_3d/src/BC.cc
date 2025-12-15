#include "../include/BC.h"

BC::BC (const GetPot& dataFile,
	    const std::string& problem,
        const std::string& section1,
        const std::string& section2) :
        M_section1 ( section1 + problem ),
        M_section2 ( section2 + problem),
	    M_nBoundaries( dataFile ( ( M_section1 + "nBoundaries" ).data (), 6) ),
	    M_BCstring( dataFile ( ( M_section1 + "bcflag" ).data (), "1") ),
	    M_BCstringF( dataFile ( ( M_section2 + "bcflag" ).data (), "1") ),
 	    M_BCNeum( dataFile ( ( M_section1 + "p_BC" ).data (), "1") ),
 	    M_BCNeumF( dataFile ( ( M_section2 + "p_BC" ).data (), "1") ),
 	    M_BCDiri( dataFile ( ( M_section1 + "v_BC" ).data (), "1") ),
 	    M_BCDiriVec( dataFile ( ( M_section1 + "bdDisp" ).data (), "1") ),
 	    M_BCNeumVec( dataFile ( ( M_section1 + "bdLoad" ).data (), "1") ),
 	    M_BCDiriVel( dataFile ( ( M_section1 + "bdVel" ).data (), "0") ),
 	    M_BCDiriF( dataFile ( ( M_section2 + "v_BC" ).data (), "1") )

{   

    M_BC.resize(M_nBoundaries,0); 

    M_parser.setString ( M_BCstring );

    
    for ( size_type i = 0; i < M_nBoundaries; ++i )
    {
        M_BC [ i ] = M_parser.evaluate ( i );

        if (M_BC[i]==0)
        {
        	M_DiriRG.push_back(i);
        }
        if (M_BC[i]==1)
        {
        	M_NeumRG.push_back(i);
        }
    	if (M_BC[i]==2)
	    {
        	M_MixedRG.push_back(i);
	    }
}

    M_BCF.resize(2,0); 
    M_parser.setString ( M_BCstringF );

    for ( size_type i = 0; i < 2;++i )
    {
        M_BCF [ i ] = M_parser.evaluate ( i );
        if (M_BCF[i]==0)
        {
        	M_DiriRGF.push_back(i);
        }
        if (M_BCF[i]==1)
        {
        	M_NeumRGF.push_back(i);
        }
    }
    M_bdNodesX.resize(M_nBoundaries+1,0);
    std::string xnodes(dataFile ( ( M_section1 + "boundaryNodesX" ).data (), "[0,1,1,0]"));
    M_parser.setString (xnodes);


    for ( size_type i = 0; i < M_nBoundaries; ++i )
    {

        M_bdNodesX [ i ] = M_parser.evaluate ( i );
     
    }
    M_bdNodesX [ M_nBoundaries ]= M_bdNodesX [ 0 ];
M_bdNodesY.resize(M_nBoundaries+1,0);
        
    std::string ynodes(dataFile ( ( M_section1 + "boundaryNodesY" ).data (), "[0,1,1,0]"));
    M_parser.setString ( ynodes );

    for ( size_type i = 0; i < M_nBoundaries; ++i )
    {
        M_bdNodesY [ i ] = M_parser.evaluate ( i );
        
    }
M_bdNodesY [ M_nBoundaries ]= M_bdNodesY [ 0 ];

   
}

scalar_type BC::BCNeum(const base_node& x, const size_type& flag, const scalar_type t)
{
    M_parser.setString ( M_BCNeum);
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    M_parser.setVariable ( "z", x [ 2 ] );
    M_parser.setVariable ( "n", flag);
    M_parser.setVariable ( "t", t);
    return M_parser.evaluate ();
}

scalar_type BC::BCNeumF(const base_node& x, const size_type& flag)
{
    M_parser.setString ( M_BCNeumF);
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    M_parser.setVariable ( "z", x [ 2 ] );
    M_parser.setVariable ( "n", flag );
    return M_parser.evaluate ();
}

scalar_type BC::BCDiri(const base_node& x, const size_type& flag)
{
    M_parser.setString ( M_BCDiri);
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    M_parser.setVariable ( "z", x [ 2 ] );
    M_parser.setVariable ( "n", flag );
    return M_parser.evaluate ();
}



bgeot::base_node BC::BCDiriVec(const base_node& x, const size_type& flag, const scalar_type t)
{
    bgeot::base_node sol(3,0);
       
    for ( size_type i = 0; i < 3; ++i )
    {
 	 M_parser.setString ( M_BCDiriVec);
   	 M_parser.setVariable ( "x", x [ 0 ] );
   	 M_parser.setVariable ( "y", x [ 1 ] );
   	 M_parser.setVariable ( "z", x [ 2 ] );
     M_parser.setVariable ( "n", flag );
	 M_parser.setVariable ( "t", t );
    	 sol[i]=M_parser.evaluate (i);
    }
    return sol;
}

bgeot::base_node BC::BCNeumVec(const base_node& x, const size_type& flag, const scalar_type t)
{
    bgeot::base_node sol(3,0);
       
    for ( size_type i = 0; i < 3; ++i )
    {
 	 M_parser.setString ( M_BCNeumVec);
   	 M_parser.setVariable ( "x", x [ 0 ] );
   	 M_parser.setVariable ( "y", x [ 1 ] );
   	 M_parser.setVariable ( "z", x [ 2 ] );
     M_parser.setVariable ( "n", flag );
	 M_parser.setVariable ( "t", t );
    	 sol[i]=M_parser.evaluate (i);
    }

    return sol;
}

bgeot::base_node BC::BCDiriVel(const base_node& x, const size_type& flag, const scalar_type t)
{
    bgeot::base_node sol(3,0);
    for ( size_type i = 0; i < 3; ++i )
    {
 	 M_parser.setString ( M_BCDiriVel);
   	 M_parser.setVariable ( "x", x [ 0 ] );
  	 M_parser.setVariable ( "y", x [ 1 ] );
  	 M_parser.setVariable ( "z", x [ 2 ] );
     M_parser.setVariable ( "n", flag );
	 M_parser.setVariable ( "t", t );
    	 sol[i]=M_parser.evaluate (i);
    }

    return sol;
}

scalar_type BC::BCDiriF(const base_node& x, const size_type& flag)
{
    M_parser.setString ( M_BCDiriF);
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    M_parser.setVariable ( "z", x [ 2 ] );
    M_parser.setVariable ( "n", flag );
    return M_parser.evaluate ();
    return 0.0;
}

std::vector<size_type> BC::getNeumBD(std::string where)
	{

		return M_NeumRG;
	}

std::vector<size_type> BC::getDiriBD(std::string where)
	{

		return M_DiriRG;
	}

std::vector<size_type> BC::getMixedBD(std::string where)
	{

		return M_MixedRG;
	}

void BC::setBoundaries(getfem::mesh* meshPtr, std::string where)
{
    getfem::mesh_region border_faces;
    getfem::outer_faces_of_mesh(*meshPtr, border_faces);
    std::vector<int> verticalBDs;
    std::vector<int> horizontalBDs;

    for (int j=0; j<M_nBoundaries;++j)
    {
			if (M_bdNodesX[j]==M_bdNodesX[j+1])
				verticalBDs.push_back(j);
			if (M_bdNodesY[j]==M_bdNodesY[j+1])
				horizontalBDs.push_back(j);			
    }


   
  	  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) 
    	  {
      		assert(i.is_face());
	
	        base_node un = meshPtr->normal_of_face_of_convex(i.cv(), i.f());
     		un /= gmm::vect_norm2(un);
            /*
            if (gmm::abs(un[0]+1)<1.0e-7 && gmm::abs(un[1])<1.0e-7 && gmm::abs(un[2])<1.0e-7)
      		{
      			meshPtr->region(0).add(i.cv(),i.f());	      	//nx=-1	
      		}
	        if (gmm::abs(un[0]-1)<1.0e-7 && gmm::abs(un[1])<1.0e-7 && gmm::abs(un[2])<1.0e-7)
      		{
      			meshPtr->region(1).add(i.cv(),i.f());	      //nx=1		
      		}
            
            if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1]+1)<1.0e-7 && gmm::abs(un[2])<1.0e-7)
      		{
      			meshPtr->region(2).add(i.cv(),i.f());	      	//ny=-1	
      		}
	        if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1]-1)<1.0e-7 && gmm::abs(un[2])<1.0e-7)
      		{
      			meshPtr->region(3).add(i.cv(),i.f());	      //ny=1		
      		}
	        
            if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1])<1.0e-7 && gmm::abs(un[2]+1)<1.0e-7)
      		{
      			meshPtr->region(4).add(i.cv(),i.f());	      	//nz=-1	
      		}
	        if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1])<1.0e-7 && gmm::abs(un[2]-1)<1.0e-7)
      		{
      			meshPtr->region(5).add(i.cv(),i.f());	      //nz=1	
      		}
            */

            //per il caso cilindrico

      		if (gmm::abs(un[2])<0.1)
      		{
      		    // Calcolo del baricentro nel piano XY
                base_node bary = gmm::mean_value(meshPtr->points_of_face_of_convex(i.cv(), i.f()));
                
                // Prodotto scalare tra normale e raggio
                double rx = bary[0] - 1.0;
                double ry = bary[1] - 1.0;
                double s = un[0]*rx + un[1]*ry;    
                if (s > 0) {
                    meshPtr->region(0).add(i.cv(), i.f());//superficie esterna
                } else {
                    meshPtr->region(1).add(i.cv(),i.f());//superficie interna
                 }	      //nz=0		
      		}
	        
            if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1])<1.0e-7 && gmm::abs(un[2]+1)<1.0e-7)
      		{
      			meshPtr->region(2).add(i.cv(),i.f());	      	//nz=-1	
      		}
	        if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1])<1.0e-7 && gmm::abs(un[2]-1)<1.0e-7)
      		{
      			meshPtr->region(3).add(i.cv(),i.f());	      //nz=1	
      		}
                    
    	  } 

}


void BC::setBoundariesFromTags(getfem::mesh* meshPtr, 
                               const std::map<std::string, size_type>& regmap,
                               std::string where)
{
    std::cout << "\n=== BC::setBoundariesFromTags ===" << std::endl;
    std::cout << "Setting boundaries from Gmsh physical tags..." << std::endl;
    
    // Get all regions that exist in the mesh
    dal::bit_vector mesh_regions = meshPtr->regions_index();
    
    std::cout << "Available regions in mesh:" << std::endl;
    for (dal::bv_visitor i(mesh_regions); !i.finished(); ++i) {
        getfem::mesh_region region = meshPtr->region(i);
        size_t count = 0;
        for (getfem::mr_visitor it(region); !it.finished(); ++it) {
            count++;
        }
        std::cout << "  Region " << i << ": " << count << " faces" << std::endl;
    }
    
    std::cout << "\nPhysical names from Gmsh:" << std::endl;
    for (const auto& pair : regmap) {
        std::cout << "  '" << pair.first << "' -> region ID " << pair.second << std::endl;
    }
    
    // Map physical names to your internal region numbering
    // This mapping can be configured in the data file

    std::map<std::string, size_type> name_to_internal_id;
    
    // Try to find standard names in regmap
    for (const auto& pair : regmap) {
        std::string name_lower = pair.first;
        // Convert to lowercase for case-insensitive matching
        std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
        
        // Map based on common naming conventions
        if (name_lower.find("outer") != std::string::npos || 
            name_lower.find("external") != std::string::npos ||
            name_lower.find("ext") != std::string::npos) {
            name_to_internal_id[pair.first] = 0;
            std::cout << "  Mapping '" << pair.first << "' to internal region 0 (outer surface)" << std::endl;
        }
        else if (name_lower.find("inner") != std::string::npos || 
                 name_lower.find("internal") != std::string::npos ||
                 name_lower.find("int") != std::string::npos) {
            name_to_internal_id[pair.first] = 1;
            std::cout << "  Mapping '" << pair.first << "' to internal region 1 (inner surface)" << std::endl;
        }
        else if (name_lower.find("bottom") != std::string::npos || 
                 name_lower.find("lower") != std::string::npos ||
                 name_lower.find("bot") != std::string::npos) {
            name_to_internal_id[pair.first] = 2;
            std::cout << "  Mapping '" << pair.first << "' to internal region 2 (bottom surface)" << std::endl;
        }
        else if (name_lower.find("top") != std::string::npos || 
                 name_lower.find("upper") != std::string::npos) {
            name_to_internal_id[pair.first] = 3;
            std::cout << "  Mapping '" << pair.first << "' to internal region 3 (top surface)" << std::endl;
        }
        else {
            std::cout << "  Warning: No mapping found for '" << pair.first << "', skipping" << std::endl;
        }
    }
    
    // Now assign the Gmsh regions to your internal regions
    std::cout << "\nAssigning Gmsh regions to internal regions..." << std::endl;
    
    for (const auto& mapping : name_to_internal_id) {
        const std::string& physical_name = mapping.first;
        size_type internal_id = mapping.second;
        
        // Find the Gmsh region ID for this physical name
        auto it = regmap.find(physical_name);
        if (it != regmap.end()) {
            size_type gmsh_region_id = it->second;
            
            // Check if this region exists in the mesh
            if (mesh_regions.is_in(gmsh_region_id)) {
                // Copy all faces from Gmsh region to internal region
                getfem::mesh_region gmsh_region = meshPtr->region(gmsh_region_id);
                
                size_t face_count = 0;
                for (getfem::mr_visitor face_it(gmsh_region); !face_it.finished(); ++face_it) {
                    meshPtr->region(internal_id).add(face_it.cv(), face_it.f());
                    face_count++;
                }
                
                std::cout << "  ✓ Assigned " << face_count << " faces from '" 
                          << physical_name << "' (Gmsh region " << gmsh_region_id 
                          << ") to internal region " << internal_id << std::endl;
            }
            else {
                std::cout << "  ✗ Warning: Gmsh region " << gmsh_region_id 
                          << " for '" << physical_name << "' not found in mesh!" << std::endl;
            }
        }
    }
    
    // Verify the assignment
    std::cout << "\nVerification - Internal regions after assignment:" << std::endl;
    for (size_type i = 0; i < M_nBoundaries; ++i) {
        size_t count = 0;
        getfem::mesh_region region = meshPtr->region(i);
        for (getfem::mr_visitor it(region); !it.finished(); ++it) {
            count++;
        }
        
        std::string bc_type = "Unknown";
        if (M_BC[i] == 0) bc_type = "Dirichlet";
        else if (M_BC[i] == 1) bc_type = "Neumann";
        else if (M_BC[i] == 2) bc_type = "Mixed";
        
        std::cout << "  Internal region " << i << ": " << count 
                  << " faces (BC type: " << bc_type << ")" << std::endl;
        
        if (count == 0) {
            std::cout << "    WARNING: No faces assigned to this region!" << std::endl;
        }
    }
    
    std::cout << "=== Boundary assignment complete ===" << std::endl;
}
