#ifndef BC_H
#define BC_H

#include "Core.h"
#include "Parser.h"

//gestione delle condizioni al contorno per i due domini (bulk e faglia). Generico rispetto al problema

class BC
{
public: 
	 BC ( const GetPot& dataFile,
	      const std::string& problem,
                          const std::string& section1 = "bulkData/",
                          const std::string& section2 = "fractureData/");
	
	scalar_type BCNeum(const base_node& x, const size_type& flag, const scalar_type t);
	scalar_type BCNeumF(const base_node& x, const size_type& flag);
	scalar_type BCDiri(const base_node& x, const size_type& flag);
	scalar_type BCDiriF(const base_node& x, const size_type& flag);
	scalar_type BCMixed(const base_node& x, const size_type& flag);  // per bordi in cui si ammette scorrimento solo parallelo
	
	bgeot::base_node BCDiriVec(const base_node& x, const size_type& flag, const scalar_type t);

	bgeot::base_node BCDiriVel(const base_node& x, const size_type& flag, const scalar_type t);
	
	bgeot::base_node BCNeumVec(const base_node& x, const size_type& flag, const scalar_type t);
	
	void setBoundaries(getfem::mesh* meshPrt,  std::string where="bulk");
	
	std::vector<size_type> getNeumBD(std::string where="bulk");
	
	std::vector<size_type> getDiriBD(std::string where="bulk");

        std::vector<size_type> getMixedBD(std::string where="bulk");
	
	
private:

    std::string M_section1;
    std::string M_section2;

    // Attributes
    std::string M_BCstring;
    std::string M_BCstringF;
    
    std::string M_BCNeum;
    std::string M_BCNeumF;
    std::string M_BCNeumVec;

    std::string M_BCDiri;
    std::string M_BCMixed;
    std::string M_BCDiriVec;
    std::string M_BCDiriVel;
    std::string M_BCDiriF;
    
    std::vector<size_type>  M_NeumRG;
    std::vector<size_type>  M_DiriRG;
    std::vector<size_type>  M_MixedRG;
    
    std::vector<size_type>  M_NeumRGF;
    std::vector<size_type>  M_DiriRGF;

    std::vector<size_type> M_BC;
    std::vector<size_type> M_BCF;
   
    LifeV::Parser M_parser;
    
    int M_nBoundaries;
    std::vector<scalar_type> M_bdNodesX;
    std::vector<scalar_type> M_bdNodesY;

};

#endif
