#ifndef BULKELASTDATA_H
#define BULKELASTDATA_H

#include "Core.h"
#include "Parser.h"

//dati meccanica relativi al bulk 

class BulkElastData
{
public: 
	 BulkElastData ( const GetPot& dataFile,
                          const std::string& section = "bulkData/",
                          const std::string& sectionDarcy = "mecc/");
	
	scalar_type Lambda(const base_node& x );
	scalar_type Mu(const base_node& x );
		
        void setLambda(std::vector<base_node> nodes);
        void setMu(std::vector<base_node> nodes);
        void setfluidP(std::vector<base_node> nodes);
         
	inline scalar_type getLambda(size_type i)
	{	
		return (*M_LambdaVector)[i];
	}
	
	inline scalar_type getMu(size_type i)
	{	
		return (*M_MuVector)[i];
	}
	
        inline std::vector<scalar_type> getLambda()
	{	
		return (*M_LambdaVector);
	}
	
	inline std::vector<scalar_type> getMu()
	{	
		return (*M_MuVector);
	}

	bgeot::base_node bulkLoad(bgeot::base_node x, scalar_type t);
	bgeot::base_node uEx(bgeot::base_node x, scalar_type t);
	bgeot::base_node uIni(bgeot::base_node x);
	scalar_type fluidP(const base_node& x );
	
	
private:

    // Attributes
    std::string M_section;
    std::string M_sectionElast;
    
    std::string M_mu;
    std::string M_lambda;
    std::string M_load;
    std::string M_uEx;
    std::string M_uIni;
    std::string M_fluidP;
    std::string M_Gstring;

    scalar_type M_rhoR;
    bgeot::base_node M_G;
        
    scalarVectorPtr_Type M_LambdaVector;
    scalarVectorPtr_Type M_MuVector;
    scalarVectorPtr_Type M_LoadX;
    scalarVectorPtr_Type M_LoadY;
    scalarVectorPtr_Type M_fluidPVector;
    
    LifeV::Parser M_parser;
};

#endif
