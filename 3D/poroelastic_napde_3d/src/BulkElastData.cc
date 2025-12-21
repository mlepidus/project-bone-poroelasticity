#include "../include/BulkElastData.h"

BulkElastData::BulkElastData ( const GetPot& dataFile,
                                 const std::string& section,
                                 const std::string& sectionElast) :
        M_section ( section ),
    	M_sectionElast ( M_section + sectionElast ),
	    M_mu( dataFile ( ( M_sectionElast + "mu" ).data (), "1." ) ),			
	    M_lambda( dataFile ( ( M_sectionElast + "lambda" ).data (), "1." ) ),	
	    M_load( dataFile ( ( M_sectionElast + "bulkLoad" ).data (), "1." ) ),
	    M_uEx( dataFile ( ( M_sectionElast + "u_exact" ).data (), "0." ) ),
	    M_uIni( dataFile ( ( M_sectionElast + "u_ini" ).data (), "0." ) ),
	    M_fluidP( dataFile ( ( M_sectionElast + "fluidP" ).data (), "0." ) ),		
	    M_Gstring( dataFile ( ( std::string("materials/") + "gravity" ).data (), "[0,0,0]" ) ),
        M_rhoR( dataFile ( ( std::string("materials/") + "rho_r" ).data (), 1000 ) )
{
     M_G.resize(3);
     M_parser.setString ( M_Gstring );


    for ( size_type i = 0; i < 3; ++i )
    {
        M_G [ i ] = M_parser.evaluate ( i );
    }
	
}

scalar_type BulkElastData::Lambda(const base_node& x )
	 {
		M_parser.setString ( M_lambda);
    		M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
		M_parser.setVariable ( "z", x [ 2 ] );
 	        return M_parser.evaluate ();
	 }

scalar_type BulkElastData::Mu(const base_node& x )
	 {
		M_parser.setString ( M_mu);
    		M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
		M_parser.setVariable ( "z", x [ 2 ] );
 	        return M_parser.evaluate ();
	 }

void BulkElastData::setLambda(std::vector<base_node> nodes)
	{

		M_LambdaVector.reset(new scalarVector_Type (nodes.size()));
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_LambdaVector)[i]=Lambda(nodes[i]);		
			}
	}
	
void BulkElastData::setMu(std::vector<base_node> nodes)
	{
		M_MuVector.reset(new scalarVector_Type (nodes.size()));
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_MuVector)[i]=Mu(nodes[i]);		
			}
	}
	
void BulkElastData::setfluidP(std::vector<base_node> nodes)
	{
		M_fluidPVector.reset(new scalarVector_Type (nodes.size()));
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_fluidPVector)[i]=fluidP(nodes[i]);		
			}
	}
	
bgeot::base_node BulkElastData::bulkLoad(bgeot::base_node x, scalar_type t)
	{
		bgeot::base_node sol(3,0);
       
	        for ( size_type i = 0; i < 3; ++i )
    		{
 	 		M_parser.setString ( M_load);
   		        M_parser.setVariable ( "x", x [ 0 ] );
    	 		M_parser.setVariable ( "y", x [ 1 ] );
			M_parser.setVariable ( "z", x [ 2 ] );
			M_parser.setVariable ( "t", t );
     	 		sol[i]=M_parser.evaluate (i);
   		 }

	        return sol;

	}

bgeot::base_node BulkElastData::uEx(bgeot::base_node x, scalar_type t)
	{
		bgeot::base_node sol(3,0);
       
	        for ( size_type i = 0; i < 3; ++i )
    		{
 	 		M_parser.setString ( M_uEx);
   		        M_parser.setVariable ( "x", x [ 0 ] );
    	 		M_parser.setVariable ( "y", x [ 1 ] );
			M_parser.setVariable ( "z", x [ 2 ] );
			M_parser.setVariable ( "t", t);
     	 		sol[i]= M_parser.evaluate (i);
   		 }

	        return sol;

	}


bgeot::base_node BulkElastData::uIni(bgeot::base_node x)
	{
		bgeot::base_node sol(3,0);
       
	        for ( size_type i = 0; i < 3; ++i )
    		{
 	 		M_parser.setString ( M_uIni);
   		        M_parser.setVariable ( "x", x [ 0 ] );
    	 		M_parser.setVariable ( "y", x [ 1 ] );
     	 		M_parser.setVariable ( "z", x [ 2 ] );
     	 		sol[i]= M_parser.evaluate (i);
   		 }

	        return sol;

	}
scalar_type BulkElastData::fluidP(const base_node& x )
	 {
		M_parser.setString ( M_fluidP);
    		M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
		M_parser.setVariable ( "z", x [ 2 ] );
 	        return M_parser.evaluate ();
	 }



