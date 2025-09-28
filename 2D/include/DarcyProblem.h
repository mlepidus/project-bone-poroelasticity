#ifndef DARCYPROBLEM_H
#define DARCYPROBLEM_H

#include "LinearSystem.h"
#include "DarcyOperators.h"
#include "UsefulFunctions.h"

//classe per il problema di darcy stazionario. contiene un puntatore al suo sistema lineare, alle geometrie, i dati e gli elementi finiti

class my_scalar_data : public getfem::nonlinear_elem_term
{
   
    bgeot::multi_index sizes_;
    std::string expression;

public:

    my_scalar_data( std::string bla);

    inline const bgeot::multi_index& sizes (bgeot::size_type ) const
    {
        return sizes_;
    }

    void compute ( getfem::fem_interpolation_context& ctx,
                           bgeot::base_tensor& t );
private: 
    mutable LifeV::Parser M_parser;
};

class my_vector_data : public getfem::nonlinear_elem_term
{
   
    bgeot::multi_index sizes_;
    std::string expression;

public:

    my_vector_data( std::string bla);

    inline const bgeot::multi_index& sizes (bgeot::size_type ) const
    {
        return sizes_;
    }

    void compute ( getfem::fem_interpolation_context& ctx,
                           bgeot::base_tensor& t );
private: 
    mutable LifeV::Parser M_parser;
};

class DarcyProblem
{
public: 
	 DarcyProblem ( const GetPot& dataFile, Bulk* bulk=NULL);
	 
	 FEM* getFEM(const std::string where="bulk", const std::string variable="Pressure");
	  
         void addToSys(LinearSystem* sys);
         
         void assembleMatrix(LinearSystem* sys);
         
         void assembleRHS(LinearSystem* sys);
         
         void solve();
         
         void extractSol(scalarVectorPtr_Type sol);
                
         void exportVtk(std::string folder="./vtk", std::string what="all");
         
         void createVtkfromFile(std::string folder, std::string what, std::string sourceFile);
         
         size_type getNDOF(std::string variable="all");

	 scalar_type computeError(std::string what="all");

	 void recontructVelocityDOF(scalarVectorPtr_Type velDOFs);
	 void recontructVelocityDOFlocal(scalarVectorPtr_Type velDOFs, size_type icv);
	 
	 inline scalarVector_Type getSol(const std::string variable="Pressure")
	 {
		if (variable=="Pressure")
			return *M_pressureSol;
		else
			return *M_velocitySol;
  	}
	inline getfem::mesh_im getIM()
	{
		return M_intMethod;
	}
	
 private:

    scalar_type error_given_function(scalarVector_Type V, getfem::mesh_fem& femP, getfem::mesh_im& im);
    scalar_type error_given_functionV(scalarVector_Type V, getfem::mesh_fem& femV, getfem::mesh_im& im);
    void L2_projection(scalarVector_Type& V, getfem::mesh_fem& femP, getfem::mesh_im& im);

    Bulk* M_Bulk;

    BC M_BC;
    FEM M_PressureFEM;

    FEM M_CoeffFEM;

    FEM M_VelocityFEM;

    FEM M_visualizationVFEM;
    
    LinearSystem* M_Sys;
    
    scalarVectorPtr_Type M_pressureSol;
    scalarVectorPtr_Type M_velocitySol;

    sparseMatrixPtr_Type A22;
    sparseMatrixPtr_Type A11;
    size_type M_nbTotDOF;
    size_type M_nbTotBulkDOF;
    
    getfem::mesh_im M_intMethod;
    
    mutable LifeV::Parser M_parser;

};


#endif

