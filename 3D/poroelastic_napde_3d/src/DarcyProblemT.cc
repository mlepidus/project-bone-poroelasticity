#include "../include/DarcyProblemT.h"

DarcyProblemT::DarcyProblemT (const GetPot& dataFile, Bulk* bulk):
			M_Bulk(bulk),
			M_BC(dataFile, "darcy/"),
			M_PressureFEM( bulk->getMesh(), dataFile, "darcy/", "Pressure", "bulkData/"),
			M_CoeffFEM( bulk->getMesh(), dataFile, "darcy/", "Coeff", "bulkData/"),
			M_VelocityFEM( bulk->getMesh(), dataFile, "darcy/", "Velocity", "bulkData/"),
			M_visualizationVFEM( bulk->getMesh(), dataFile, "darcy/", "VelocityVis", "bulkData/"),
			M_intMethod(*(bulk->getMesh()) ),
			M_Sys()

{  
   M_nbTotDOF=M_PressureFEM.nb_dof()+M_VelocityFEM.nb_dof();
  
   std::string intMethod(dataFile ( std::string("bulkData/darcy/integrationMethod" ).data (), "IM_TETRAHEDRON(2)" ) );
   M_intMethod.set_integration_method(bulk->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod) );
   
   M_Bulk->getDarcyData()->setKxx(M_CoeffFEM.getDOFpoints());
   M_Bulk->getDarcyData()->setKyy(M_CoeffFEM.getDOFpoints());
   M_Bulk->getDarcyData()->setKxy(M_CoeffFEM.getDOFpoints());
   M_Bulk->getDarcyData()->setKzz(M_CoeffFEM.getDOFpoints());
   M_Bulk->getDarcyData()->setKxz(M_CoeffFEM.getDOFpoints());
   M_Bulk->getDarcyData()->setKyz(M_CoeffFEM.getDOFpoints());

      if (M_Bulk->hasExternalMesh()) {
       std::cout << "Using Gmsh physical tags for boundary conditions..." << std::endl;
       M_BC.setBoundariesFromTags(M_Bulk->getMesh(), M_Bulk->getRegionMap(), "bulk");
   } else {
       std::cout << "Using geometric detection for boundary conditions..." << std::endl;
       M_BC.setBoundaries(M_Bulk->getMesh(), "bulk");
   }
   
}

FEM* DarcyProblemT::getFEM(const std::string where, const std::string variable)
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

size_type DarcyProblemT::getNDOF(std::string variable)
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
		return M_VelocityFEM.nb_dof()+M_PressureFEM.nb_dof();
	}
	
}
void DarcyProblemT::addToSys(LinearSystem* sys)
{
	M_Sys=sys;
        M_Sys->addToMatrix(M_nbTotDOF);
   
   
}

void DarcyProblemT::initialize()
{
	M_pressureSolOld.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
	M_pressureSol.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
	
	gmm::clear(*M_pressureSol);
	
	M_velocitySol.reset(new scalarVector_Type (M_VelocityFEM.nb_dof()));
	
	gmm::clear(*M_velocitySol);
	
	for (size_type i=0; i<M_PressureFEM.nb_dof("base");++i)
	{
		base_node nodo(M_PressureFEM.point_of_basic_dof(i));
		(*M_pressureSol)[i]=M_Bulk->getDarcyData()->pIni(nodo);
		(*M_pressureSolOld)[i]=(*M_pressureSol)[i];
	}


}

void DarcyProblemT::updateSol()
{	
	gmm::clear(*M_pressureSolOld);
	gmm::copy(*M_pressureSol, *M_pressureSolOld);
}

void DarcyProblemT::assembleMatrix(LinearSystem* sys, std::string where)
{
		sparseMatrixPtr_Type A11;
		A11.reset(new sparseMatrix_Type (M_VelocityFEM.nb_dof(),M_VelocityFEM.nb_dof()));
		gmm::clear(*A11);
		sparseMatrixPtr_Type A12;
		A12.reset(new sparseMatrix_Type (M_VelocityFEM.nb_dof(),M_PressureFEM.nb_dof()));
		gmm::clear(*A12);
		if (A22==NULL)
		{
		A22.reset(new sparseMatrix_Type (M_PressureFEM.nb_dof(),M_PressureFEM.nb_dof()));
		gmm::clear(*A22);
		}
		massHdiv( A11, M_Bulk, M_VelocityFEM, M_CoeffFEM, M_intMethod);
 		massL2( A22, M_Bulk, M_PressureFEM, M_CoeffFEM, M_intMethod, M_dt);
				
		divHdiv( A12, M_Bulk, M_VelocityFEM, M_PressureFEM, M_intMethod);

		essentialWNitsche( A11, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod);
	
		M_Sys->addSubMatrix(A11, 0,0);
			

		M_Sys->addSubMatrix(A12, 0, gmm::mat_ncols(*A11));
		M_Sys->addSubMatrix(A12, gmm::mat_nrows(*A11),0, -1.0, true);
		if (M_Bulk->getDarcyData()->M()>0)
		M_Sys->addSubMatrix(A22, M_VelocityFEM.nb_dof(),M_VelocityFEM.nb_dof(),1./M_Bulk->getDarcyData()->M());

		M_Sys->saveMatrix("mat_darcy.mm");

}

void DarcyProblemT::assembleRHS(LinearSystem* sys, std::string where)
{
	
			scalarVectorPtr_Type  source;
			source.reset(new scalarVector_Type (M_PressureFEM.getFEM()->nb_dof()));
			gmm::clear(*source);

			scalarSource(source, M_Bulk, M_PressureFEM, M_CoeffFEM, M_intMethod, M_timeLoop->time() );
		
			M_Sys->addSubVector(source,M_VelocityFEM.nb_dof());
		
			scalarVectorPtr_Type  BCvec;
			BCvec.reset(new scalarVector_Type (M_VelocityFEM.getFEM()->nb_dof()));
			gmm::clear(*BCvec);
			essentialWNitscheRHS( BCvec, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod);
			M_Sys->addSubVector(BCvec, 0);		
		
			gmm::clear(*BCvec);
		
			naturalRHS( BCvec, M_Bulk, &M_BC, M_VelocityFEM, M_CoeffFEM, M_intMethod, M_timeLoop->time() );
			//std::cout << "bc fatte"<<std::endl;

			M_Sys->addSubVector(BCvec, 0);	

			gmm::clear(*BCvec);

                        vectorSource( BCvec, M_Bulk, M_VelocityFEM, M_CoeffFEM, M_intMethod);
		
			M_Sys->addSubVector(BCvec, 0);
			if (M_Bulk->getDarcyData()->M()>0)
			M_Sys->multAddToRHS(A22,M_pressureSolOld, 0, M_VelocityFEM.nb_dof(),1./M_Bulk->getDarcyData()->M());



}

void DarcyProblemT::solve()
{
       	M_Sys->solve();
       	M_pressureSol.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
	    gmm::clear(*M_pressureSol);
       	M_Sys->extractSubVector(M_pressureSol, M_VelocityFEM.nb_dof(), "sol");
 
       	gmm::clear(*M_velocitySol);
      	M_Sys->extractSubVector(M_velocitySol, 0, "sol");

	if (M_timeLoop->time()==0)
	{
	    M_pressureSolIni.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
		M_Sys->extractSubVector(M_pressureSolIni, M_VelocityFEM.nb_dof(), "sol");
	}

}

void DarcyProblemT::extractSol(scalarVectorPtr_Type sol)
{
	M_pressureSol.reset(new scalarVector_Type (M_PressureFEM.nb_dof()));
	gmm::clear(*M_pressureSol);
	gmm::copy ( gmm::sub_vector (*sol, gmm::sub_interval ( M_VelocityFEM.nb_dof(),  M_PressureFEM.nb_dof()) ), *M_pressureSol);

    M_velocitySol.reset(new scalarVector_Type (M_VelocityFEM.nb_dof()));
	gmm::clear(*M_velocitySol);
    gmm::copy ( gmm::sub_vector (*sol, gmm::sub_interval ( 0, M_VelocityFEM.nb_dof()) ), *M_velocitySol);
       	
 
}

scalar_type DarcyProblemT::computeError(std::string what, scalar_type time)
{
	std::cout << "compute error"<<std::endl;
	scalar_type error;
	scalarVector_Type loc_err(M_PressureFEM.getFEM()->nb_dof());
	if (what=="Pressure" || what=="all")
	{
		for (int i=0;i<loc_err.size();++i)
		{
			bgeot::base_node where=M_PressureFEM.getFEM()->point_of_basic_dof(i);
			loc_err[i]=gmm::abs(((*M_pressureSol)[i]) - M_Bulk->getDarcyData()->pEx(where,time) );
		}
	}	
	error=L2Norm ( A22,loc_err, M_Bulk, *(M_PressureFEM.getFEM()), M_intMethod, -1);
	error=error*M_timeLoop->dt();
	std::cout <<"at time   "<<time<< "   error pressure   "<<error<<std::endl; 
	return error;
}

void DarcyProblemT::exportVtk(std::string folder, std::string what, int frame)
{
	std::cout << "export Darcy solution"<<std::endl;
	std::string frameNum;
	
	if (frame<0)
	{
		frameNum="";
	}
	else
	{
		if (frame<10)
		{
			frameNum.append("00"+LifeV::number2string(frame));
		}
		else
		{
			if (frame<100)
			{
				frameNum.append("0"+LifeV::number2string(frame));
			}
			else
			{
				frameNum.append(LifeV::number2string(frame));
			}
		}
	}		

	
 	
	if (what=="Pressure" || what=="all")
	{
		
		getfem::vtk_export exp(folder + "Pressure_bulk"+ frameNum + ".vtk" );
		getfem::vtk_export expE(folder + "Pressure_ex"+ frameNum + ".vtk" );
		getfem::vtk_export expD(folder + "dP"+ frameNum + ".vtk" );
  		
		
			exp.exporting( *(M_PressureFEM.getFEM()));
			exp.write_mesh();
			expE.exporting( *(M_PressureFEM.getFEM()));
			expE.write_mesh();
			expD.exporting( *(M_PressureFEM.getFEM()));
			expD.write_mesh();
			std::vector<scalar_type> diff(M_PressureFEM.getFEM()->nb_dof());
			std::vector<scalar_type> pe(M_PressureFEM.getFEM()->nb_dof());
			for (int i=0; i<diff.size(); ++i)
			{
				bgeot::base_node where=M_PressureFEM.getFEM()->point_of_basic_dof(i);
				pe[i]=M_Bulk->getDarcyData()->pEx(where,M_timeLoop->time() );
			}
			if (M_pressureSolIni!=NULL)
			{
			for (int i=0; i<diff.size(); ++i)
			{
				diff[i]=(*M_pressureSol)[i] - (*M_pressureSolIni)[i] ;
							}
			}
			if (true)//(frame>0)
			{
			exp.write_point_data( *(M_PressureFEM.getFEM()), *M_pressureSol, "p");
			expE.write_point_data( *(M_PressureFEM.getFEM()), pe, "pex");
			expD.write_point_data( *(M_PressureFEM.getFEM()), diff, "pIni");
			}
			else
			{
			exp.write_point_data( *(M_PressureFEM.getFEM()), *M_pressureSolOld, "p");
			}
		
		
		
	}
	
	if (what=="Velocity" || what=="all")
	{
		getfem::vtk_export  exp(folder + "Velocity_bulk"+ frameNum + ".vtk" );
		
	
			std::vector<scalar_type> velInt(M_visualizationVFEM.getFEM()->nb_dof(),0.0);
			getfem::interpolation(*(M_VelocityFEM.getFEM()), *(M_visualizationVFEM.getFEM()),*M_velocitySol, velInt );
  			exp.exporting(*(M_visualizationVFEM.getFEM()));
			exp.write_mesh();
			exp.write_point_data( *(M_visualizationVFEM.getFEM()), velInt, "v");
		
		
		
		
	}
	
}

