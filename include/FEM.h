#ifndef FEM_H
#define FEM_H

#include "Core.h"


//questa classe vuole essere un wrapper per la classe FEM di getfem e renderla un po' più friendly

class FEM
{
public: 
	 FEM (const getfem::mesh* mesh,
	      const GetPot& dataFile,
	      const std::string& problem,
	      const std::string& variable,
              const std::string& section = "bulkData/"
              );
              
 	FEM (const getfem::mesh* mesh,
	         std::string femType,
   		 size_type spaceDim
              );

         size_type nb_dof(std::string which="all");    //numero gradi di libertà: posso chiedere quelli base, estesi o totali

	 inline std::string type()
	 {
		return M_femType;
	 }  
         
         inline getfem::mesh_fem* getFEM()   //puntatore al mesh_fem
         {
         	return &M_FEM;
         }
         inline std::vector<base_node> getDOFpoints()  //punti corrispondenti ai dof
         {	
         	return M_DOFpoints;
         }
         
        
         inline std::vector<size_type> getExt()  //restituisce vettore con gli identificatori dei dof raddoppiati
         {
         	return M_extended;
         }
         
         inline std::vector<scalar_type> getExtSign()  //il dof è dalla parte + o - del level set?
         {
         	return M_extSign;
         }
         
         inline size_type getExt(size_type i)  //l'i-esimo dof esteso a quale base corrisponde?
         {
         	return M_extended[i];
         }
         
	 inline scalar_type getDOFSign(size_type i)
	 {
		return M_DOFSign[i];
	 }

         inline scalar_type getExtSign(size_type i)
         {
         	return M_extSign[i];
         } 
         
         inline base_node point_of_basic_dof(size_type i)
         {
         	return M_FEM.point_of_basic_dof(i);
         }
          
private:

    std::string M_section;
    std::string M_femType;
    size_type M_SpaceDim;
    getfem::mesh_fem M_FEM;
    const getfem::mesh* M_meshPtr;
    std::vector<base_node> M_DOFpoints;
 //   mutable LifeV::Parser M_parser;
    
    std::vector<size_type> M_extended;
    std::vector<scalar_type> M_extSign;
    std::vector<scalar_type> M_DOFSign;
};

#endif
