#include "../include/CoupledProblem.h"

CoupledProblem::CoupledProblem (const GetPot& dataFile, Bulk* bulk,  TimeLoop* time ):
			M_Bulk(bulk),
			M_intMethod(*(bulk->getMesh()) ),
		        M_time(time),
			M_Sys()

{  
   std::string intMethod(dataFile ( std::string("bulkData/darcy/integrationMethod" ).data (), "IM_TRIANGLE(6)" ) );
   M_intMethod.set_integration_method(bulk->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod) );
}


void CoupledProblem::addToSys(LinearSystem* sys)
{ 
	M_Sys=sys;
	   std::cout << "n dof darcy"<< M_DarcyPB->getNDOF()<< "   n dof elast"<<M_ElastPB->getNDOF()<<std::endl;
        M_Sys->addToMatrix(M_ElastPB->getNDOF()+M_DarcyPB->getNDOF());
       
        M_ElastPB->addToSys(&M_elastSys);
        M_DarcyPB->addToSys(&M_darcySys);
           
}

void CoupledProblem::assembleMatrix(LinearSystem* sys)
{
	M_ElastPB->assembleMatrix(&M_elastSys,"bulk");

	
	M_DarcyPB->assembleMatrix(&M_darcySys,"bulk");
	
	
	M_PressureStress.reset(new sparseMatrix_Type (M_ElastPB->getFEM("bulk","Disp")->nb_dof(), M_DarcyPB->getFEM("bulk","Pressure")->nb_dof()));
	
	matrixFluidP( M_PressureStress, M_Bulk, *(M_ElastPB->getFEM("bulk","Disp")), (*M_DarcyPB->getFEM("bulk","Pressure")), M_intMethod);	
}

void CoupledProblem::assembleRHS(LinearSystem* sys)
{
	M_ElastPB->assembleRHS(&M_elastSys,"bulk");

	M_DarcyPB->assembleRHS(&M_darcySys,"bulk");

	//M_DarcyPB->solve();
}

void CoupledProblem::addSubSystems()
{
	M_Sys->addSubSystem(&M_elastSys,0,0);

	M_Sys->addSubSystem(&M_darcySys,M_ElastPB->getNDOF(),M_ElastPB->getNDOF());

	M_Sys->addSubMatrix(M_PressureStress, 0,M_ElastPB->getNDOF()+ M_DarcyPB->getNDOF("Velocity"));
	
	M_Sys->addSubMatrix(M_PressureStress, M_ElastPB->getNDOF()+ M_DarcyPB->getNDOF("Velocity"),0, -1.0/(*M_time).dt(), true);
	
	M_Sys->multAddToRHS(M_ElastPB->getOldSol(), M_ElastPB->getNDOF()+ M_DarcyPB->getNDOF("Velocity"),0 ,M_DarcyPB->getNDOF("Pressure"),  M_ElastPB->getNDOF("Disp") );
	
	
}

void CoupledProblem::addSubSystemsRHS()
{
	M_Sys->addSubSystemRHS(&M_elastSys,0,0);

	M_Sys->addSubSystemRHS(&M_darcySys,M_ElastPB->getNDOF(),M_ElastPB->getNDOF());

	M_Sys->multAddToRHS(M_ElastPB->getOldSol(), M_ElastPB->getNDOF()+ M_DarcyPB->getNDOF("Velocity"),0 ,M_DarcyPB->getNDOF("Pressure"),  M_ElastPB->getNDOF("Disp") );
	
	
}
void CoupledProblem::clearSubSystems()
{
	gmm::clear(*(M_darcySys.getRHS()));
	gmm::clear((*M_elastSys.getRHS()));
	gmm::clear((*M_darcySys.getMatrix()));
	gmm::clear((*M_elastSys.getMatrix()));

}
void CoupledProblem::clearSubSystemsRHS()
{
	gmm::clear(*(M_darcySys.getRHS()));
	gmm::clear((*M_elastSys.getRHS()));
}

void CoupledProblem::solve()
{
	M_Sys->solve();
       	M_solEl.reset(new scalarVector_Type (M_ElastPB->getNDOF()));
	gmm::clear(*M_solEl);
       	M_Sys->extractSubVector(M_solEl, 0, "sol");
   	M_solDa.reset(new scalarVector_Type (M_DarcyPB->getNDOF()));
	gmm::clear(*M_solDa);
       	M_Sys->extractSubVector(M_solDa, M_ElastPB->getNDOF(), "sol");
       	
       	M_ElastPB->extractSol(M_solEl);
       	M_DarcyPB->extractSol(M_solDa);

}

void CoupledProblem::updateSol()
{
	M_ElastPB->updateSol();
       	M_DarcyPB->updateSol();
 

}

void CoupledProblem::enforceStrongBC(bool firstTime)
{
    if (firstTime)
    {
    	
	M_ElastPB->enforceStrongBC(true);
    }
    else
   {
	M_ElastPB->enforceStrongBC(false);

    }

}


bgeot::base_node CoupledProblem::computeError(scalar_type t)
{
	bgeot::base_node error(0.0,0.0);
	scalar_type ep=M_DarcyPB->computeError("all",t);
	error[0]=ep;
	return error;
	
}

void CoupledProblem::exportVtk(std::string folder, std::string what, int frame)
{
	if (what=="all")
	{
		M_DarcyPB->exportVtk(folder, "all",frame);

		M_ElastPB->exportVtk(folder, "all",frame);
	}
	else
	{
		if (what=="elast")
		{
			M_ElastPB->exportVtk(folder, "all",frame);
		}
		else
		{
			M_DarcyPB->exportVtk(folder, "all",frame);
		}
	}

}
