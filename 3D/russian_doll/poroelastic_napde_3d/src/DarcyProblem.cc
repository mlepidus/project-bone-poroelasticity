#include "../include/DarcyProblem.h"

my_scalar_data:: my_scalar_data (  std::string bla) :
     expression(bla)
{
    sizes_.resize(1);
    sizes_ [ 0 ] = 1;
 }
void my_scalar_data::compute ( getfem::fem_interpolation_context& ctx,
                                      bgeot::base_tensor& t )
{
    size_type cv = ctx.convex_num(); 
    bgeot::base_node where=ctx.xreal();
    M_parser.setString ( expression);
    M_parser.setVariable ( "x", where [ 0 ] );
    M_parser.setVariable ( "y", where [ 1 ] );
    M_parser.setVariable ( "t", 0 );
    t [ 0 ] = M_parser.evaluate();
}

my_vector_data:: my_vector_data (  std::string bla) :
     expression(bla)
{
    sizes_.resize(1);
    sizes_ [ 0 ] = 2;
 }

void my_vector_data::compute ( getfem::fem_interpolation_context& ctx,
                                      bgeot::base_tensor& t )
{
    size_type cv = ctx.convex_num(); 
    bgeot::base_node where=ctx.xreal();
    for (size_type i=0;i<2;++i)
    {
    M_parser.setString ( expression);
    M_parser.setVariable ( "x", where [ 0 ] );
    M_parser.setVariable ( "y", where [ 1 ] );
    M_parser.setVariable ( "t", 0 );
    t [ i ] = M_parser.evaluate(i);
    }
}


DarcyProblem::DarcyProblem (const GetPot& dataFile, Bulk* bulk):
			M_Bulk(bulk),
			M_BC(dataFile, "darcy/"),
			M_PressureFEM( bulk->getMesh(), dataFile, "darcy/", "Pressure", "bulkData/"),
			M_CoeffFEM( bulk->getMesh(), dataFile, "darcy/", "Coeff", "bulkData/"),
			M_VelocityFEM( bulk->getMesh(), dataFile, "darcy/", "Velocity", "bulkData/"),
			M_visualizationVFEM( bulk->getMesh(), dataFile, "darcy/", "VelocityVis", "bulkData/"),
         	M_intMethod(*(bulk->getMesh()) ),
			M_Sys()

{  
  
   M_nbTotBulkDOF=M_PressureFEM.nb_dof()+M_VelocityFEM.nb_dof();
   M_nbTotDOF=M_PressureFEM.nb_dof()+M_VelocityFEM.nb_dof();
   
   std::string intMethod(dataFile ( std::string("bulkData/darcy/integrationMethod" ).data (), "IM_TETRAHEDRON(2)" ) );
   M_intMethod.set_integration_method(bulk->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod) );
      
   M_Bulk->getDarcyData()->setKxx(M_CoeffFEM.getDOFpoints());
   M_Bulk->getDarcyData()->setKyy(M_CoeffFEM.getDOFpoints());
   M_Bulk->getDarcyData()->setKxy(M_CoeffFEM.getDOFpoints());
   M_BC.setBoundaries(M_Bulk->getMesh());  
}


FEM* DarcyProblem::getFEM(const std::string where, const std::string variable)
 {
 	if (variable=="Pressure")
 	{
 		
			return &M_PressureFEM;
	}	
	
 	if (variable=="Velocity")
 	{
 	
			return &M_VelocityFEM;
	}
}

size_type DarcyProblem::getNDOF(std::string variable)
{
	if (variable=="Pressure")
	{
		return M_PressureFEM.nb_dof();
	}

	if (variable=="Velocity")
	{
		return M_VelocityFEM.nb_dof();
	}
     
	if (variable=="all")
	{
		return M_PressureFEM.nb_dof()+M_VelocityFEM.nb_dof();
	}
	
}


void DarcyProblem::addToSys(LinearSystem* sys)
{
	M_Sys=sys;
        M_Sys->addToMatrix(M_nbTotDOF);
}


void DarcyProblem::assembleMatrix(LinearSystem* sys)
{
			sparseMatrixPtr_Type A11;
			A11.reset(new sparseMatrix_Type (M_VelocityFEM.nb_dof(),M_VelocityFEM.nb_dof()));
			gmm::clear(*A11);
			sparseMatrixPtr_Type A12;
			A12.reset(new sparseMatrix_Type (M_VelocityFEM.nb_dof(),M_PressureFEM.nb_dof()));
			gmm::clear(*A12);
			
			massHdiv( A11, M_Bulk, M_VelocityFEM, M_CoeffFEM, M_intMethod);
				
			divHdiv( A12, M_Bulk, M_VelocityFEM, M_PressureFEM, M_intMethod);
			
			gmm::MatrixMarket_IO::write("B.mm" , *A12); 

			essentialWNitsche( A11, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod);
	
			M_Sys->addSubMatrix(A11, 0,0);
			M_Sys->addSubMatrix(A12, 0, gmm::mat_ncols(*A11));
			M_Sys->addSubMatrix(A12, gmm::mat_nrows(*A11),0, -1.0, true);

			A22.reset(new sparseMatrix_Type (M_PressureFEM.nb_dof(),M_PressureFEM.nb_dof()));
			gmm::clear(*A22);
			
			massL2( A22, M_Bulk, M_PressureFEM, M_CoeffFEM, M_intMethod, 1.);
			

}
void DarcyProblem::assembleRHS(LinearSystem* sys)
{
	
			scalarVectorPtr_Type  source;
			source.reset(new scalarVector_Type (M_PressureFEM.getFEM()->nb_dof()));
			gmm::clear(*source);
		
			scalarSource(source, M_Bulk, M_PressureFEM, M_CoeffFEM, M_intMethod,0);
				
			M_Sys->addSubVector(source,M_VelocityFEM.nb_dof());
		
			scalarVectorPtr_Type  BCvec;
			BCvec.reset(new scalarVector_Type (M_VelocityFEM.getFEM()->nb_dof()));
			gmm::clear(*BCvec);

			essentialWNitscheRHS( BCvec, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod);

			M_Sys->addSubVector(BCvec, 0);		
		
			gmm::clear(*BCvec);
		
			naturalRHS( BCvec, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod,0);

			M_Sys->addSubVector(BCvec, 0);	

			gmm::clear(*BCvec);

			vectorSource( BCvec, M_Bulk, M_VelocityFEM, M_CoeffFEM, M_intMethod);
			
			M_Sys->addSubVector(BCvec, 0);


}

void DarcyProblem::solve()
{
       	M_Sys->solve();
       	M_pressureSol.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
	gmm::clear(*M_pressureSol);
       	M_Sys->extractSubVector(M_pressureSol, M_VelocityFEM.nb_dof(), "sol");
 
       	M_velocitySol.reset(new scalarVector_Type (M_VelocityFEM.nb_dof()));
	gmm::clear(*M_velocitySol);
      	M_Sys->extractSubVector(M_velocitySol, 0, "sol");
      	
}

void DarcyProblem::extractSol(scalarVectorPtr_Type sol)
{
	M_pressureSol.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
	gmm::clear(*M_pressureSol);
	gmm::copy ( gmm::sub_vector (*sol, gmm::sub_interval ( M_VelocityFEM.nb_dof(),  M_PressureFEM.nb_dof()) ), *M_pressureSol);

       	M_velocitySol.reset(new scalarVector_Type (M_VelocityFEM.nb_dof()));
	gmm::clear(*M_velocitySol);
        gmm::copy ( gmm::sub_vector (*sol, gmm::sub_interval ( 0, M_VelocityFEM.nb_dof()) ), *M_velocitySol);
    
}


scalar_type DarcyProblem::computeError(std::string what)
{
   	
	std::cout << "compute error"<<std::endl;
	scalar_type error1;
	scalar_type error2;
	scalar_type error3;
	
	if (what=="Pressure" || what=="all")
	{
		scalarVector_Type loc_err(M_PressureFEM.getFEM()->nb_dof());
		for (int i=0;i<loc_err.size();++i)
		{
			bgeot::base_node where=M_PressureFEM.getFEM()->point_of_basic_dof(i);
			loc_err[i]=gmm::abs(((*M_pressureSol)[i]) -M_Bulk->getDarcyData()->pEx(where,0) );
			
		}
		error1=L2Norm ( A22,loc_err, M_Bulk, *(M_PressureFEM.getFEM()), M_intMethod, -1);

		error2=error_given_function(*M_pressureSol,*(M_PressureFEM.getFEM()), M_intMethod);
	}	

	if (what=="Velocity" || what=="all")
	{
		error3=error_given_functionV(*M_velocitySol,*(M_VelocityFEM.getFEM()), M_intMethod);
	}	
		
	
	std::cout <<"at time   "<<time<< "   error pressure   "<<error1<<std::endl; 
	std::cout <<"at time   "<<time<< "   error pressure (metodo cattivo)   "<<error2<<std::endl; 
	std::cout <<"at time   "<<time<< "   error velocity   "<<error3<<std::endl; 
	return error1;
}

scalar_type DarcyProblem::error_given_function(scalarVector_Type V, getfem::mesh_fem& femP, getfem::mesh_im& im)
{
	getfem::generic_assembly assem;
	std::vector<scalar_type> v(femP.nb_dof());

	my_scalar_data nterm( M_Bulk->getDarcyData()->getPexpr()); 
	assem.set("u=data(#1);"
          "V(#1)+=u(i).u(j).comp(Base(#1).Base(#1).Base(#1))(:,i,j)-2*u(j).comp(Base(#1).Base(#1).NonLin(#2))(:,j,1)+comp(Base(#1).NonLin(#2).NonLin(#2))(:,1,1)");
        assem.push_mi(im);       
        assem.push_mf(femP);
        assem.push_mf(femP);
	assem.push_nonlinear_term(&nterm);
	assem.push_data(V);
	assem.push_vec(v);
	assem.assembly(-1);
	scalar_type somma=0;
	for (size_type i=0; i<femP.nb_dof();++i)
	{
		somma+=v[i];
	}
	return pow(somma,0.5);
}

scalar_type DarcyProblem::error_given_functionV(scalarVector_Type V, getfem::mesh_fem& femV, getfem::mesh_im& im)
{
	getfem::generic_assembly assem;
	std::vector<scalar_type> v(1);
	my_vector_data nterm(M_Bulk->getDarcyData()->getUexpr());
	assem.set("u=data(#1);"
          "V()+=u(i).u(j).comp(vBase(#1).vBase(#1))(i,k,j,k)+comp(NonLin(#2).NonLin(#2))(i,i)-u(i).comp(vBase(#1).NonLin(#2))(i,k,k)-u(i).comp(vBase(#1).NonLin(#2))(i,k,k)");
        assem.push_mi(im);       
        assem.push_mf(femV);
        assem.push_mf(femV);
	assem.push_nonlinear_term(&nterm);
	assem.push_data(V);
	assem.push_vec(v);
	assem.assembly(-1);
	return pow(v[0],0.5);
}

void DarcyProblem::L2_projection(scalarVector_Type& V, getfem::mesh_fem& femP, getfem::mesh_im& im)
{
	getfem::generic_assembly assem;
	std::vector<scalar_type> aa(V.size());

	for (size_type i=0; i<femP.nb_dof();++i)
	{
		size_type ii=femP.first_convex_of_basic_dof(i);
		aa[ii]=1.0/M_Bulk->getMesh()->convex_area_estimate(i);
	}

	my_scalar_data nterm( M_Bulk->getDarcyData()->getPexpr());
//  
	assem.set("a=data(#1);"
          "V(#1)+=a(i).comp(Base(#1).Base(#1).NonLin(#2))(:,i,1)");
        assem.push_mi(im);       
        assem.push_mf(femP);
        assem.push_mf(femP);
	assem.push_nonlinear_term(&nterm);
	assem.push_data(aa);
	assem.push_vec(V);
	assem.assembly(-1);
}


void DarcyProblem::recontructVelocityDOF(scalarVectorPtr_Type velDOFs)
{
	scalarVectorPtr_Type  projSol;
	projSol.reset(new scalarVector_Type (M_VelocityFEM.getFEM()->nb_dof()));
	gmm::clear(*projSol);
	scalarVectorPtr_Type  provv;
	provv.reset(new scalarVector_Type (M_visualizationVFEM.getFEM()->nb_dof()));
	gmm::clear(*provv);

	for (int i=0; i<M_visualizationVFEM.getFEM()->nb_dof(); ++++i)
	{
		bgeot::base_node where=M_visualizationVFEM.getFEM()->point_of_basic_dof(i);
		bgeot::base_node uex=M_Bulk->getDarcyData()->uEx(where,0);
		(*provv)[i]=uex[0];
		(*provv)[i+1]=uex[1];
	}

	vectorSource (projSol, M_Bulk,  M_VelocityFEM, M_visualizationVFEM, provv, M_intMethod);
	sparseMatrixPtr_Type mass;
	mass.reset(new sparseMatrix_Type (M_VelocityFEM.nb_dof(),M_VelocityFEM.nb_dof()));
	gmm::clear(*mass);

	massL2( mass,  M_Bulk,   M_VelocityFEM, M_CoeffFEM, M_intMethod, 1);
	scalar_type rcond;
  	gmm::clear(*velDOFs);
	SuperLU_solve(*mass, *velDOFs, *projSol, rcond);
	
}
void DarcyProblem::recontructVelocityDOFlocal(scalarVectorPtr_Type velDOFs, size_type icv)
{
	scalarVectorPtr_Type  projSol;
	projSol.reset(new scalarVector_Type (M_VelocityFEM.getFEM()->nb_dof()));
	gmm::clear(*projSol);
	scalarVectorPtr_Type  provv;
	provv.reset(new scalarVector_Type (M_visualizationVFEM.getFEM()->nb_dof()));
	gmm::clear(*provv);



	for (int i=0;i<M_visualizationVFEM.getFEM()->ind_basic_dof_of_element(icv).size(); ++++i)
	{	
		size_type quale=M_visualizationVFEM.getFEM()->ind_basic_dof_of_element(icv)[i];
		bgeot::base_node where=M_visualizationVFEM.getFEM()->point_of_basic_dof(quale);
		bgeot::base_node uex=M_Bulk->getDarcyData()->uEx(where,0);
		(*provv)[quale]=uex[0];
		(*provv)[quale+1]=uex[1];
	}

	vectorSource (projSol, M_Bulk,  M_VelocityFEM, M_visualizationVFEM, provv, M_intMethod);
	sparseMatrixPtr_Type mass;
	mass.reset(new sparseMatrix_Type (M_VelocityFEM.nb_dof(),M_VelocityFEM.nb_dof()));
	gmm::clear(*mass);

	massL2Local( mass,  M_Bulk,   M_VelocityFEM, M_CoeffFEM, M_intMethod, 1,icv);
	scalar_type rcond;
  	gmm::clear(*velDOFs);
	SuperLU_solve(*mass, *velDOFs, *projSol, rcond);
	
}

void DarcyProblem::exportVtk(std::string folder, std::string what)
{
	std::cout << "export"<<std::endl;
	
	
 	
	if (what=="Pressure" || what=="all")
	{
		getfem::vtk_export exp(folder + "Pressure_bulk.vtk");
  		
		
			exp.exporting( *(M_PressureFEM.getFEM()));
			exp.write_mesh();
			exp.write_point_data( *(M_PressureFEM.getFEM()), *M_pressureSol, "p");
		

	}

	if (what=="Velocity" || what=="all")
	{
		getfem::vtk_export exp(folder + "Velocity_bulk.vtk");
		
		
			std::vector<scalar_type> velInt(M_visualizationVFEM.getFEM()->nb_dof(),0.0);
			getfem::interpolation(*(M_VelocityFEM.getFEM()), *(M_visualizationVFEM.getFEM()),*M_velocitySol, velInt );
  			exp.exporting(*(M_visualizationVFEM.getFEM()));
			exp.write_mesh();
			exp.write_point_data( *(M_visualizationVFEM.getFEM()), velInt, "v");
				
	}
	
}

void DarcyProblem::createVtkfromFile(std::string folder, std::string what, std::string sourceFile)
{
	std::cout << "create vtk from file"<<std::endl;
	
	std::ifstream myfile; myfile.open(folder+sourceFile);

        
 	
	if (what=="Pressure" )
	{
	    std::vector<scalar_type> v(M_PressureFEM.getFEM()->nb_dof(),0.0);
	    for (size_type i=0; i<v.size();++i)
	    {
	       myfile>>v[i];

	    }
		getfem::vtk_export exp(folder + "p_from_file.vtk");
  		
		
			exp.exporting( *(M_PressureFEM.getFEM()));
			exp.write_mesh();
			exp.write_point_data( *(M_PressureFEM.getFEM()), v, "p");
		

	}

	if (what=="Velocity" )
	{
	    std::vector<scalar_type> v(M_VelocityFEM.getFEM()->nb_dof(),0.0);
	    for (size_type i=0; i<v.size();++i)
	    {
	       myfile>>v[i];

	    }
		getfem::vtk_export exp(folder + "v_from_file.vtk");
		
		exp.exporting( *(M_VelocityFEM.getFEM()));
			exp.write_mesh();
			exp.write_point_data( *(M_VelocityFEM.getFEM()), v, "v");
				
	}
	
}
