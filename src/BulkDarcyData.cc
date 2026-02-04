#include "../include/BulkDarcyData.h"

BulkDarcyData::BulkDarcyData ( const GetPot& dataFile,
                                 const std::string& section,
                                 const std::string& sectionDarcy) :
            M_section ( section ),
            M_sectionDarcy ( M_section + sectionDarcy ),
	    M_Kxx( dataFile ( ( M_sectionDarcy + "Kxx" ).data (), "1." ) ),
 	    M_Kyy( dataFile ( ( M_sectionDarcy + "Kyy" ).data (), "1." ) ),
	    M_Kxy( dataFile ( ( M_sectionDarcy + "Kxy" ).data (), "1." ) ),
		M_Kzz( dataFile ( ( M_sectionDarcy + "Kzz" ).data (), "1." ) ),
 	    M_Kxz( dataFile ( ( M_sectionDarcy + "Kxz" ).data (), "1." ) ),
	    M_Kyz( dataFile ( ( M_sectionDarcy + "Kyz" ).data (), "1." ) ),
	    M_source( dataFile ( ( M_sectionDarcy + "source" ).data (), "1." ) ),
	    M_pIni( dataFile ( ( M_sectionDarcy + "pIni" ).data (), "0." ) ),
		M_pBC( dataFile ( ( M_sectionDarcy + "p_BC" ).data (), "0." ) ),
	    M_pEx( dataFile ( ( M_sectionDarcy + "p_exact" ).data (), "0." ) ),
	    M_uEx( dataFile ( ( M_sectionDarcy + "u_exact" ).data (), "[0,0,0]" ) ),
	    M_Gstring( dataFile( (std::string("materials/") + "gravity").data(), "[0,0,0]" ) ),
	    M_rhoF( dataFile( (std::string("materials/") + "rho_f").data(), 1000 ) ),
	    M_muF( dataFile( (std::string("materials/") + "mu_f").data(), 1.0e-3 ) ),
 	    M_biotM( dataFile( (M_sectionDarcy + "biotM" ).data(), 1.0 ) ),
		M_leakage( dataFile( (M_sectionDarcy + "leakage" ).data(), 1.0 ) ),
		M_biotAlpha( dataFile( (M_sectionDarcy + "biotAlpha" ).data(), 1.0 ) )
{
	M_G.resize(3);
    M_parser.setString ( M_Gstring );
	for (size_type i=0; i<3; ++i)
	{
	    M_G [ i ] = M_parser.evaluate ( i );
	}

}


scalar_type BulkDarcyData::Kxx(const base_node& x )
	 {	
		size_type dim = x.size();
		M_parser.setString ( M_Kxx);
    	M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
 	    if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
 	    return M_parser.evaluate ()/M_muF;
	 }

scalar_type BulkDarcyData::Kxy(const base_node& x )
	 {	
		size_type dim = x.size();
		M_parser.setString ( M_Kxy);
    	M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
		if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
 	    return M_parser.evaluate ()/M_muF;
	 }

scalar_type BulkDarcyData::Kyy(const base_node& x )
	 {	
		size_type dim = x.size(); 
		M_parser.setString ( M_Kyy);
    	M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
 	    if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
 	    return M_parser.evaluate ()/M_muF;
	 }
scalar_type BulkDarcyData::Kzz(const base_node& x )
	 {	
		size_type dim = x.size(); 
		M_parser.setString ( M_Kzz);
   		M_parser.setVariable ( "x", x [ 0 ] );
    	M_parser.setVariable ( "y", x [ 1 ] );
	    if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
 	    return M_parser.evaluate ()/M_muF;
	 }
scalar_type BulkDarcyData::Kxz(const base_node& x )
	 {	
		size_type dim = x.size(); 
		M_parser.setString ( M_Kxz);
    	M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
		if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
 	    return M_parser.evaluate ()/M_muF;
	 }
scalar_type BulkDarcyData::Kyz(const base_node& x )
	 {	
		size_type dim = x.size(); 
		M_parser.setString ( M_Kyz);
    	M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
		if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
 	    return M_parser.evaluate ()/M_muF;
	 }
	 
scalar_type BulkDarcyData::source(const base_node& x, const scalar_type t=0)
	 {	
		size_type dim = x.size(); 
		M_parser.setString ( M_source);
  		M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
 	    if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
		M_parser.setVariable ( "t", t );
 	    return M_parser.evaluate ();
	 }
	 
scalar_type BulkDarcyData::pIni(const base_node& x )
	 {	
		size_type dim = x.size(); 
		M_parser.setString ( M_pIni);
   		M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
 	    if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
 	    return M_parser.evaluate ();
	 }

scalar_type BulkDarcyData::p_BC(const base_node& x, const scalar_type t)
{  
   size_type dim = x.size();
   M_parser.setString ( M_pBC );
   M_parser.setVariable ( "x", x [ 0 ] );
   M_parser.setVariable ( "y", x [ 1 ] );
   if (dim >= 3) 
   		M_parser.setVariable("z", x[2]);
   else         
   		M_parser.setVariable("z", 0.0);
   M_parser.setVariable ( "t", t );
   return M_parser.evaluate ();
}

scalar_type BulkDarcyData::pEx(const base_node& x , const scalar_type t=0)
	 {	
		static bool first = true;
		if (first) {
			std::cout << "--- [BulkData] Parsing pEx string: " << M_pEx << " ---" << std::endl;
			first = false;
		}
		size_type dim = x.size(); 
		M_parser.setString ( M_pEx);
    	M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
		if (dim >= 3)
        	M_parser.setVariable("z", x[2]);
    	else
        	M_parser.setVariable("z", 0.0);
		M_parser.setVariable ( "t", t );
 	    return M_parser.evaluate ();
	 }


bgeot::base_node BulkDarcyData::uEx(const base_node& x, const scalar_type t=0)
	 {	
		size_type dim = x.size();
		bgeot::base_node u(dim,0);
		for (size_type i=0;i<dim;++i)
		{
			M_parser.setString ( M_uEx);
			M_parser.setVariable ( "x", x [ 0 ] );
			M_parser.setVariable ( "y", x [ 1 ] );
			if (dim >= 3)
            	M_parser.setVariable ( "z", x[2] );
        	else
            	M_parser.setVariable ( "z", 0.0 );
			M_parser.setVariable ( "t", t );
			u[i]=M_parser.evaluate (i);
		}
 	        return u;
	 }

void BulkDarcyData::setKxx(std::vector<base_node> nodes)
	{
		M_KxxVector.reset(new scalarVector_Type (nodes.size()));
		
		
		for (size_type i=0;i<nodes.size();++i)
		{
		      bgeot::base_node x(nodes[i]);

		      (*M_KxxVector)[i]=Kxx(nodes[i]);
		}
		
	}
	
void BulkDarcyData::setKyy(std::vector<base_node> nodes)
	{
		M_KyyVector.reset(new scalarVector_Type (nodes.size()));
		
		for (size_type i=0;i<nodes.size();++i)
			{
			bgeot::base_node x(nodes[i]);

			(*M_KyyVector)[i]=Kyy(nodes[i]);		
			}
		
	}
	
void BulkDarcyData::setKxy(std::vector<base_node> nodes)
	{
		M_KxyVector.reset(new scalarVector_Type (nodes.size()));
		
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_KxyVector)[i]=Kxy(nodes[i]);		
			}
		
	}

void BulkDarcyData::setKzz(std::vector<base_node> nodes)
	{
		M_KzzVector.reset(new scalarVector_Type (nodes.size()));
		
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_KzzVector)[i]=Kzz(nodes[i]);		
			}
		
	}

void BulkDarcyData::setKxz(std::vector<base_node> nodes)
	{
		M_KxzVector.reset(new scalarVector_Type (nodes.size()));
		
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_KxzVector)[i]=Kxz(nodes[i]);		
			}
		
	}

void BulkDarcyData::setKyz(std::vector<base_node> nodes)
	{
		M_KyzVector.reset(new scalarVector_Type (nodes.size()));
		
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_KyzVector)[i]=Kyz(nodes[i]);		
			}
		
	}
