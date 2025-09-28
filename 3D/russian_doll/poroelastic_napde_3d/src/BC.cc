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
 	    M_BCDiriF( dataFile ( ( M_section2 + "v_BC" ).data (), "1") ),
        M_hasOuterBoundaryInterpolationD(false),
        M_hasOuterBoundaryInterpolationP(false)
{   

    M_BC.resize(M_nBoundaries,0); 

    M_parser.setString ( M_BCstring );

    
    for ( size_type i = 0; i < M_nBoundaries; ++i )
    {
        M_BC [ i ] = M_parser.evaluate ( i );

        if (M_BC[i]==0)
        {
        	M_DiriRG.push_back(i);
            std::cout << "Dirichlet BC on boundary " << i << std::endl;
        }
        if (M_BC[i]==1)
        {
        	M_NeumRG.push_back(i);
            std::cout << "Neumann BC on boundary " << i << std::endl;
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
// Check if this is the outer boundary (region 0) and we have an interpolation function
    if (flag == 0 && M_hasOuterBoundaryInterpolationP) {
        // Use the polynomial interpolation for the outer boundary
        // Assuming the boundary condition depends on the z-coordinate (vertical direction)
        return M_outerBoundaryInterpolationP(x[2]);
    } else {
        // For all other boundaries, use the original parser-based approach
    M_parser.setString ( M_BCNeum);
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    M_parser.setVariable ( "z", x [ 2 ] );
    M_parser.setVariable ( "n", flag);
    M_parser.setVariable ( "t", t);
    return M_parser.evaluate ();
    }
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
    // Check if this is the outer boundary (region 0) and we have an interpolation function
    if (flag == 0 && M_hasOuterBoundaryInterpolationD) {
        // Use the polynomial interpolation for each component
        // Here we assume the interpolation applies to all components equally
        // Modify if different components need different handling
        sol = M_outerBoundaryInterpolationD(x[2]);
    } else {
        // Original parser-based approach for other boundaries
        for (size_type i = 0; i < 3; ++i)
        {
            M_parser.setString(M_BCDiriVec);
            M_parser.setVariable("x", x[0]);
            M_parser.setVariable("y", x[1]);
            M_parser.setVariable("z", x[2]);
            M_parser.setVariable("n", flag);
            M_parser.setVariable("t", t);
            sol[i] = M_parser.evaluate(i);
        }
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
    
                double rx = bary[0] - 1.0;
                double ry = bary[1] - 1.0;
                double s = un[0]*rx + un[1]*ry;
                // Prodotto scalare tra normale e raggio
    
                if (s > 0) {
                    meshPtr->region(0).add(i.cv(), i.f());//superficie esterna
                    //std::cout << "Outer boundary face added to region 0" << std::endl;
                } else {
                    meshPtr->region(1).add(i.cv(),i.f());//superficie interna
                   // std::cout << "Inner boundary face added to region 1" << std::endl;
                 }	      //nz=0		
      		}
	        
            if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1])<1.0e-7 && gmm::abs(un[2]+1)<1.0e-7)
      		{
      			meshPtr->region(2).add(i.cv(),i.f());	      	//nz=-1	
              //  std::cout << "bottom boundary face added to region 2" << std::endl;
      		}
	        if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1])<1.0e-7 && gmm::abs(un[2]-1)<1.0e-7)
      		{
      			meshPtr->region(3).add(i.cv(),i.f());	      //nz=1	
               // std::cout << "top boundary face added to region 3" << std::endl;
      		}
            
 
    	  } 
    
}

void BC::setOuterBoundaryInterpolationP(const std::string& filename, int degree) {
    std::vector<PointData> data = read_csv(filename);
    std::function<double(double)> x=poly_fit(data, degree);
    if (!std::isnan(x(0.0))) {
        M_outerBoundaryInterpolationP = x;
        std::cout << "Outer boundary pressure interpolated" << std::endl;
        M_hasOuterBoundaryInterpolationP = true;
    }
    else{
       std::cout << "Error in polynomial fitting: NaN encountered. Discarded interpolation." << std::endl;
    }
}

void BC::setOuterBoundaryInterpolationD(const std::string& filename, int degree) {
    std::vector<PointData> datax, datay, dataz;
    read_csv_vector(filename, datax, datay, dataz);

    auto interpX = poly_fit(datax, degree);
    auto interpY = poly_fit(datay, degree);
    auto interpZ = poly_fit(dataz, degree);

    if (!std::isnan(interpZ(0.0))&&!std::isnan(interpY(0.0))&&!std::isnan(interpX(0.0))) {
M_outerBoundaryInterpolationD = [interpX, interpY, interpZ](double z) {
    // Validate input
    if (std::isnan(z) || std::isinf(z)) {
        std::cerr << "Warning: Invalid z value in interpolation: " << z << std::endl;
        bgeot::base_node result(3, 0.0);
        return result;
    }
    
    bgeot::base_node result(3);
    
    // Evaluate each component with error handling
    try {
        result[0] = interpX(z);
        result[1] = interpY(z);
        result[2] = interpZ(z);
    } catch (const std::exception& e) {
        std::cerr << "Error in interpolation evaluation: " << e.what() << std::endl;
        result[0] = 0.0;
        result[1] = 0.0;
        result[2] = 0.0;
    }
    
    // Check for NaN in results
    if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2])) {
        std::cerr << "Warning: NaN detected in interpolation result" << std::endl;
        result[0] = 0.0;
        result[1] = 0.0;
        result[2] = 0.0;
    }
    
    return result;
};
    } else {
        std::cout << "Error in polynomial fitting for displacement: NaN encountered. Discarded interpolation." << std::endl;
    }
    std::cout << "Outer boundary displacement interpolated" << std::endl;
    M_hasOuterBoundaryInterpolationD = true;
}