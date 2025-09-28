#include "../include/FixedStress.h"


FixedStress::FixedStress (const GetPot& dataFile, Bulk* bulk, TimeLoop* time ):
			M_Bulk(bulk),
			M_intMethod(*(bulk->getMesh()) ),
		        M_time(time),
			M_betaString( dataFile ( std::string("FS/beta" ).data (), "CL" )),
			M_contactString( dataFile ( std::string("FS/contactAlgorithm" ).data (), "UZ" )),
			M_maxit( dataFile ( std::string("FS/maxit" ).data (), 100 )),
			M_toll( dataFile ( std::string("FS/toll" ).data (), 1.0e-6)),
			M_ContactToll( dataFile ( std::string("FS/ContactToll" ).data (), 1.0e-2)),
			M_ContactMaxIt( dataFile ( std::string("FS/ContactMaxIt" ).data (), 100)),
			M_isCoupled( dataFile ( std::string("time/coupled" ).data (), true))

{  
   std::string intMethod(dataFile ( std::string("bulkData/darcy/integrationMethod" ).data (), "IM_TRIANGLE(6)" ) );
   M_intMethod.set_integration_method(bulk->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod) );

   M_beta.resize(M_Bulk->getElastData()->getLambda().size());

   if (M_betaString=="CL")  //diverse scelte per il parametro beta del fixed stress
   {
	for (int i=0;i<M_beta.size();++i)
	{
		M_beta[i]=1.0/(M_Bulk->getElastData()->getLambda(i)+M_Bulk->getElastData()->getMu(i));
	}
   }
   if (M_betaString=="OPT")
   {
	for (int i=0;i<M_beta.size();++i)
	{
		M_beta[i]=0.5/(M_Bulk->getElastData()->getLambda(i)+M_Bulk->getElastData()->getMu(i));
	}
   }
   if (M_betaString=="LAMBDA")
   {
	for (int i=0;i<M_beta.size();++i)
	{
		M_beta[i]=0.5/(M_Bulk->getElastData()->getLambda(i));
	}
   }

   if (!M_isCoupled)
	std::cout << "decoupled"<<std::endl;  //se decoupled farà in pratica una sola iterazione

}

//azzera tutte le soluzioni e assembla le matrici (quelle non cambiano)

void FixedStress::initialize()
{

    M_solElOld.reset(new scalarVector_Type (M_ElastPB->getNDOF()));
    M_solDaOld.reset(new scalarVector_Type (M_DarcyPB->getNDOF()));
    M_solElPre.reset(new scalarVector_Type (M_ElastPB->getNDOF()));
    M_solDaPre.reset(new scalarVector_Type (M_DarcyPB->getNDOF()));
    M_solEl.reset(new scalarVector_Type (M_ElastPB->getNDOF()));
    M_solDa.reset(new scalarVector_Type (M_DarcyPB->getNDOF()));

    M_darcySys=new LinearSystem();
    M_elastSys=new LinearSystem();
 
    M_DarcyPB->addToSys(M_darcySys);

    M_ElastPB->addToSys(M_elastSys);

    assembleMatrices();

}

// chiama gli assembler dei due problemi nelle varie regioni

void FixedStress::assembleMatrices()
{

	//----------darcy----------------------

	M_DarcyPB->assembleMatrix(M_darcySys,"bulk");

	
       //-----------elast-----------------------

	M_ElastPB->assembleMatrix(M_elastSys,"bulk");

	//M_ElastPB->getSys()->computeInverse();

	if (M_isCoupled)
	{
     	    //------------termini aggiuntivi in darcy nella matrice-------------

     	    size_type shift=M_DarcyPB->getNDOF("Velocity");
	    size_type ssize=M_DarcyPB->getNDOF("Pressure");

            M_massMatrix.reset(new sparseMatrix_Type (ssize,ssize));

	    gmm::copy(gmm::scaled(gmm::sub_matrix(*(M_darcySys->getMatrix()), gmm::sub_interval(shift,ssize),  
			gmm::sub_interval(shift,ssize) ), M_Bulk->getDarcyData()->M()),  *M_massMatrix);
            M_darcySys->addSubMatrix(M_massMatrix, M_DarcyPB->getNDOF("Velocity"), M_DarcyPB->getNDOF("Velocity"), M_beta[0]);

            shift=M_DarcyPB->getNDOF("Velocity")+M_DarcyPB->getNDOF("Pressure")+M_DarcyPB->getNDOF("VelocityF");
            ssize=M_DarcyPB->getNDOF("PressureF");

            M_PressureStress.reset(new sparseMatrix_Type (M_ElastPB->getFEM("bulk","Disp")->nb_dof(), M_DarcyPB->getFEM("bulk","Pressure")->nb_dof()));

	    // matrici per il calcolo del contributo di pressione nel problema meccanico (non so se leuso ancora)

	    matrixFluidP( M_PressureStress, M_Bulk, *(M_ElastPB->getFEM("bulk","Disp")), (*M_DarcyPB->getFEM("bulk","Pressure")), M_intMethod);
	}
}

//risolvo sottoproblema fluido modificato dal fixed stress

void FixedStress::solveFluidStep()
{
	M_darcySys->cleanRHS();
	//assemblo rhs

	M_DarcyPB->assembleRHS(M_darcySys,"bulk");

	if (M_isCoupled)
	{
		scalarVectorPtr_Type diffU;
        	diffU.reset(new scalarVector_Type (M_ElastPB->getNDOF("Disp")));
       		gmm::clear(*diffU);

		for (int i=0;i<diffU->size();++i)
		{
			(*diffU)[i]= -(*M_solElPre)[i]+(*M_solElOld)[i];

		}
    		//aggiungo il termine di derivata temporale della divergenza di u
        	M_darcySys->multAddToRHS(M_PressureStress, diffU, 0, M_DarcyPB->getNDOF("Velocity"), -1.0/(*M_time).dt(), true);
        	M_darcySys->multAddToRHS(M_massMatrix, *M_solDaPre, M_DarcyPB->getNDOF("Velocity"), M_DarcyPB->getNDOF("Velocity"), M_beta[0] );
	}
	M_DarcyPB->solve();  //risolvo

	gmm::copy(*(M_darcySys->getSol()) , *M_solDa);
}

//risolvo sottoproblema meccanico ma SENZA gli algoritmi per il contatto.

void FixedStress::solveMechStep()
{
	M_elastSys->cleanRHS();  //riassemblo RHS
	M_ElastPB->assembleRHS(M_elastSys,"bulk");
	
	if (M_isCoupled)	
	{
		M_elastSys->multAddToRHS(M_PressureStress, *M_solDa,M_DarcyPB->getNDOF("Velocity"),0,-1.0);  //aggiungo il contributo di pressione NOTA a dx
	}
        M_ElastPB->solve();  //risolvo

	gmm::copy(*(M_elastSys->getSol()) , *M_solEl);

}

//questo contiene il vero e proprio ciclo

void FixedStress::solve()
{

	gmm::copy(*M_solElOld,*M_solElPre);
	gmm::copy(*M_solDaOld,*M_solDaPre);

	scalar_type err=M_toll+1;
	size_type cont=0;

        int mmaxit=M_maxit;

	if (!M_isCoupled)  //se non è accoppiato fa una sola iterazione
		mmaxit=1;



	while ((err>M_toll) && cont<mmaxit)
	{
	
	      solveFluidStep();
	      solveMechStep();  //non basta più questo solver, devo chiamare quelli con il problema di contatto
	     
	    gmm::copy(*(M_elastSys->getSol()) , *M_solEl);
	    err=stoppingCriterion();

	    gmm::copy(*M_solEl,*M_solElPre);
	    gmm::copy(*M_solDa,*M_solDaPre);
	    cont=cont+1;
	    std::cout <<err<<"   "<< cont<<std::endl;
	
	}	// fine while loop

	std::cout <<"convergenza"<<err<<"   "<< cont<<std::endl;
	FILE * pFile;

	//esporto i dati relativi alla convergenza
 	if (M_time->time()==0)
	{
  		pFile = fopen ( "FS_iterations.txt" , "w" );
	}	
	else
	{
  		pFile = fopen ( "FS_iterations.txt" , "a" );
	}
	fprintf ( pFile,"%e%s%d\n",M_time->time(),"  ", int(cont));
 	fclose (pFile);

}

//criterio di arresto
scalar_type FixedStress::stoppingCriterion()
{

	std::vector<scalar_type> errP(M_DarcyPB->getNDOF("Pressure"),0);

	for (int i=0;i<errP.size();++i)	
	{
		errP[i]=gmm::abs((*M_solDa)[i+M_DarcyPB->getNDOF("Velocity")]-(*M_solDaPre)[i+M_DarcyPB->getNDOF("Velocity")]);
	}

        std::vector<scalar_type> errD(M_ElastPB->getNDOF("Disp"),0);
	for (int i=0;i<errD.size();++i)	
	{
		errD[i]=gmm::abs((*M_solEl)[i]-(*M_solElPre)[i]);
	}
	std::vector<scalar_type>::const_iterator it,it2,it3,it4;
        
	it = std::max_element(errP.begin(), errP.end());
	it2 = std::max_element(errD.begin(), errD.end());
 	it3 = std::max_element((*M_solDa).begin()+M_DarcyPB->getNDOF("Velocity"), (*M_solDa).end());
	it4 = std::max_element((*M_solEl).begin(), (*M_solEl).end());
	if (*it3!=0 && *it4!=0)
	{
	return *it/gmm::abs(*it3)+*it2/gmm::abs(*it4);
	}
	else
	{
	return *it + *it2;
	}
}


void FixedStress::updateSol()
{
	gmm::copy(*M_solEl,*M_solElOld);
 	gmm::copy(*M_solDa,*M_solDaOld);

	M_ElastPB->updateSol();
       	M_DarcyPB->updateSol();
}


void FixedStress::exportVtk(std::string folder, std::string what, int frame)
{
	if (what=="all")
	{
		M_DarcyPB->exportVtk(folder, "all",frame);
		M_ElastPB->exportVtk(folder, "Displacement",frame);
	}
	else
	{
		if (what=="elast")
		{
			M_ElastPB->exportVtk(folder, "Displacement",frame);
		}
		else
		{
			M_DarcyPB->exportVtk(folder, "all",frame);
		}
	}
	
}
