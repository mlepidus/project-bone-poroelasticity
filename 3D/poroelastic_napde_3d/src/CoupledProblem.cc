#include "../include/CoupledProblem.h"

CoupledProblem::CoupledProblem (const GetPot& dataFile, Bulk* bulk,  TimeLoop* time,const std::string& section):
			M_Bulk(bulk),
			M_time(time),
			M_section(section),
    		M_Sys(nullptr),
      		M_ElastPB(nullptr),
     		M_DarcyPB(nullptr),
			M_intMethod(*(bulk->getMesh()) ),
			step(0)

{  
   // Rileva la dimensione della mesh caricata
   size_type meshDim = bulk->getMesh()->dim();
   
   // Se la mesh è 2D usa triangoli, se 3D usa tetraedri
   std::string defaultMethod = (meshDim == 2) ? "IM_TRIANGLE(2)" : "IM_TETRAHEDRON(2)";

   // Cerca nel file, se non trova usa il default "intelligente"
   std::string intMethod(dataFile ( std::string(M_section + "darcy/integrationMethod" ).data (), defaultMethod.c_str() ) );
   
   M_intMethod.set_integration_method(bulk->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod) );
}


void CoupledProblem::addToSys(LinearSystem* sys)
{ 
	M_Sys=sys;
	   std::cout << "n dof darcy:"<< M_DarcyPB->getNDOF()<< "   n dof elast:"<<M_ElastPB->getNDOF()<<std::endl;
        M_Sys->addToMatrix(M_ElastPB->getNDOF()+M_DarcyPB->getNDOF());
       
        M_ElastPB->addToSys(&M_elastSys);
        M_DarcyPB->addToSys(&M_darcySys);
           
}

void CoupledProblem::assembleMatrix()
{
	M_ElastPB->assembleMatrix();

	
	M_DarcyPB->assembleMatrix();
	
	
	M_PressureStress.reset(new sparseMatrix_Type (M_ElastPB->getFEM()->nb_dof(), M_DarcyPB->getFEM("Pressure")->nb_dof()));
	
	matrixFluidP( M_PressureStress, *(M_ElastPB->getFEM()), (*M_DarcyPB->getFEM("Pressure")), M_intMethod);	
}

void CoupledProblem::assembleRHS()
{
	M_ElastPB->assembleRHS();

	M_DarcyPB->assembleRHS();

	//M_DarcyPB->solve();
}

void CoupledProblem::addSubSystems()
{
	M_Sys->addSubSystem(&M_elastSys,0,0);

	M_Sys->addSubSystem(&M_darcySys,M_ElastPB->getNDOF(),M_ElastPB->getNDOF());
	
	scalar_type alpha = M_Bulk->getDarcyData()->getBiotAlpha();
    scalar_type dt = M_time->dt();
	//1,3
	M_Sys->addSubMatrix(M_PressureStress, 0, M_ElastPB->getNDOF()+ M_DarcyPB->getNDOF("Velocity"), -alpha);
	//3,1
	M_Sys->addSubMatrix(M_PressureStress, M_ElastPB->getNDOF()+ M_DarcyPB->getNDOF("Velocity"),0, alpha/dt, true);
	
	M_Sys->multAddToRHS(M_ElastPB->getOldSol(), M_ElastPB->getNDOF()+ M_DarcyPB->getNDOF("Velocity"),0 ,M_DarcyPB->getNDOF("Pressure"),  M_ElastPB->getNDOF("Disp") );
	
	
}

void CoupledProblem::addSubSystemsRHS()
{
	M_Sys->addSubSystemRHS(&M_elastSys,0);

	M_Sys->addSubSystemRHS(&M_darcySys,M_ElastPB->getNDOF());

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

	step++;
	// Esporta i dati del time step corrente
	exportHistory(*M_solEl, "displacement_history.csv", step);
	exportHistory(M_DarcyPB->getPressureSolution(), "pressure_history.csv", step);
	exportHistory(M_DarcyPB->getVelocitySolution(), "velocity_history.csv", step);

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
	scalar_type eu=M_ElastPB->computeError("Displacement",t);
	error[0]=ep;
	error[1]=eu;
	return error;
	
}

void CoupledProblem::exportVtk(std::string folder, std::string what, int frame)
{
	if (what=="all")
	{
		M_DarcyPB->exportVtk(folder, "all",frame);

		M_ElastPB->exportVtk(folder,frame);
	}
	else
	{
		if (what=="elast")
		{
			M_ElastPB->exportVtk(folder,frame);
		}
		else
		{
			M_DarcyPB->exportVtk(folder, "all",frame);
		}
	}

}

void CoupledProblem::exportHistory(const scalarVector_Type& timestepData, const std::string& filename, size_t step)
{
    // Apri il file in modalità append
    std::ofstream outFile(filename, std::ios::app);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // Scrivi i dati del time step nel file
    outFile << "Time Step " << step << ",";
    for (size_t i = 0; i < timestepData.size(); ++i)
    {
        outFile << timestepData[i];
        if (i < timestepData.size() - 1)
            outFile << ",";
    }
    outFile << "\n";

    outFile.close();
    std::cout << "Data for time step " << step << " appended to " << filename << std::endl;
}