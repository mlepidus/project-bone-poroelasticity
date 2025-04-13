#ifndef BULKDARCYDATA_H
#define BULKDARCYDATA_H

#include "Core.h"
#include "Parser.h"

//dati problema fluido nel bulk


class BulkDarcyData
{
public: 
	 BulkDarcyData ( const GetPot& dataFile,
                          const std::string& section = "bulkData/",
                          const std::string& sectionDarcy = "darcy/");
	
	scalar_type Kxx(const base_node& x );
	scalar_type Kyy(const base_node& x );
        scalar_type Kxy(const base_node& x );
        
        void setKxx(std::vector<base_node> nodes);
        void setKxy(std::vector<base_node> nodes);
        void setKyy(std::vector<base_node> nodes);

        bgeot::base_node spaceDistrK(const base_node& x );

	inline scalar_type getKxx(size_type i)
	{	
		return (*M_KxxVector)[i];
	}
	inline scalar_type getKxy(size_type i)
	{
		return (*M_KxyVector)[i];
	}
	inline scalar_type getKyy(size_type i)
	{
		return (*M_KyyVector)[i];
	}

        inline scalar_type rhoF()
	{
		return M_rhoF;
	}

        inline bgeot::base_node gravity()
	{
		return M_G;
	}

	scalar_type source(const base_node& x, const scalar_type t);
        scalar_type pIni(const base_node& x );
        scalar_type pEx(const base_node& x, const scalar_type t );
	bgeot::base_node uEx(const base_node& x, const scalar_type t);
	inline scalar_type M()
	{
		return M_biotM;
	}

	void getLayers();

	inline std::string getPexpr()
	{
		return M_pEx;
	}
        inline std::string getUexpr()
	{
		return M_uEx;
	}
        
private:

    // Attributes
    std::string M_section;
    std::string M_sectionDarcy;
    std::string M_Kxx;
    std::string M_Kyy;
    std::string M_Kxy;
    std::string M_source;
    std::string M_pIni;
    std::string M_pEx;
    std::string M_uEx;
    std::string M_Gstring;
    
    scalarVectorPtr_Type M_KxxVector;
    scalarVectorPtr_Type M_KyyVector;
    scalarVectorPtr_Type M_KxyVector;

    bgeot::base_node M_G;
    scalar_type M_rhoF;
    scalar_type M_muF;

    scalar_type M_biotM;
    
    scalarVectorPtr_Type M_;

    bool M_KfromLayers;
    int M_nLayers;
    std::vector<scalar_type> M_horizons;
    std::string M_horizonString;

    std::string M_KLayerXString;
    std::string M_KLayerYString;

    std::vector<scalar_type> M_KLayersX;
    std::vector<scalar_type> M_KLayersY;
    
    mutable LifeV::Parser M_parser;
};

#endif
