#include "../include/ElastProblem.h"

ElastProblem::ElastProblem (const GetPot& dataFile, Bulk* bulk ):
			M_Bulk(bulk),
			M_BC(dataFile, "mecc/"),
			M_DispFEM( bulk->getMesh(), dataFile, "mecc/", "Displacement", "bulkData/"),
			M_DispScalarFEM(bulk->getMesh(),dataFile ( std::string("bulkData/mecc/FEMTypeDisplacement").data (), "FEM_PK(3,1)" ), 1 ),
			M_StressFEM( bulk->getMesh(), dataFile, "mecc/", "Stress", "bulkData/"),
			M_StressScalarFEM( bulk->getMesh(), dataFile, "mecc/", "StressS", "bulkData/"),
			M_CoeffFEM( bulk->getMesh(), dataFile, "mecc/", "Coeff", "bulkData/"),
			M_intMethod(*(bulk->getMesh()) ),
			M_staticMu (dataFile ( (std::string("bulkData/mecc/static_mu" )).data (), 0.4) ), 
			M_IsNitSym(dataFile ( (std::string("bulkData/nitsche/isSymmetric" )).data (), true )  ),
			M_IsNitCons(dataFile ( (std::string("bulkData/nitsche/isConsistent" )).data (), true )  ),
			M_Sys()
{  


   M_nbTotDOF=M_DispFEM.nb_dof();// if I use Lagrange multipliers  I need more + M_StressFEM.nb_dof("extended"); 
   std::string intMethod(dataFile ( std::string("bulkData/mecc/integrationMethod" ).data (), "IM_TETRAHEDRON(2)" ) );
   M_intMethod.set_integration_method(bulk->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod) );


   M_Bulk->getElastData()->setLambda(M_CoeffFEM.getDOFpoints());  //valuta i coefficienti negli elementi della mesh

   M_Bulk->getElastData()->setMu(M_CoeffFEM.getDOFpoints());

   M_Bulk->getElastData()->setfluidP(M_CoeffFEM.getDOFpoints());

      if (M_Bulk->hasExternalMesh()) {
       std::cout << "Using Gmsh physical tags for boundary conditions..." << std::endl;
       M_BC.setBoundariesFromTags(M_Bulk->getMesh(), M_Bulk->getRegionMap());
   } else {
       std::cout << "Using geometric detection for boundary conditions..." << std::endl;
       M_BC.setBoundaries(M_Bulk->getMesh());
   }
}



//restituisce un puntatore ai FEM dl displacement. In teoria dovrei estenderla se avessi Lagr. mult

FEM* ElastProblem::getFEM(const std::string where, const std::string variable)
{
 		return &M_DispFEM;
}

// dimensiona il sistema lineare per contenere il problema
void ElastProblem::addToSys(LinearSystem* sys)
{
	M_Sys=sys;
        M_Sys->addToMatrix(M_nbTotDOF);

}

//conto dei dof
size_type ElastProblem::getNDOF(std::string variable)
{
	if (variable=="Disp")
	{
		return M_DispFEM.nb_dof();
	}
	
	if (variable=="Lagrange")
	{
		return M_StressFEM.nb_dof("extended");
	}
	if (variable=="all")
	{
		return M_DispFEM.nb_dof();//+M_StressFEM.nb_dof("extended");  //obsoleto perché non uso più Lagr. Mult.
	}
}

//inizializza la soluzione
void ElastProblem::initialize()
{
         	M_DispSol.reset(new scalarVector_Type (M_DispFEM.nb_dof()));
		M_DispSolOld.reset(new scalarVector_Type (M_DispFEM.nb_dof()));
         	gmm::clear(*M_DispSol);
		int counter=0;
		for (size_type i=0; i<M_DispFEM.nb_dof("base");++++i)
		{
			base_node nodo(M_DispFEM.point_of_basic_dof(i));
			(*M_DispSol)[counter]=M_Bulk->getElastData()->uIni(nodo)[0];
			(*M_DispSol)[counter+1]=M_Bulk->getElastData()->uIni(nodo)[1];
		}
}

// assembla la matrice
void ElastProblem::assembleMatrix(LinearSystem* sys, std::string where)
{
	
			sparseMatrixPtr_Type A;
			A.reset(new sparseMatrix_Type (M_DispFEM.nb_dof(),M_DispFEM.nb_dof()));
			gmm::clear(*A);

			if (!M_massMatrix)
			{
				M_massMatrix.reset(new sparseMatrix_Type (M_StressScalarFEM.nb_dof(),M_StressScalarFEM.nb_dof()));
				gmm::clear(*M_massMatrix);
			}

			if (!M_massMatrixVector)
			{
				M_massMatrixVector.reset(new sparseMatrix_Type (M_StressFEM.nb_dof(),M_StressFEM.nb_dof()));
				gmm::clear(*M_massMatrix);
			}

			if (!M_dispMassMatrix)
			{
				M_dispMassMatrix.reset(new sparseMatrix_Type (M_DispFEM.nb_dof(), M_DispFEM.nb_dof()));
				gmm::clear(*M_dispMassMatrix);
				massMatrix(M_dispMassMatrix, M_Bulk, M_DispFEM, M_intMethod);
			}
		

			stiffElast( A, M_Bulk, M_DispFEM, M_CoeffFEM, M_intMethod);
			massMatrix( M_massMatrixVector, M_Bulk, M_StressFEM, M_intMethod);
			massMatrix( M_massMatrix, M_Bulk, M_StressScalarFEM, M_intMethod);

			//essentialWNitscheVec( A,M_Bulk,  &M_BC,  M_DispFEM, M_CoeffFEM, M_intMethod);
				
			M_Sys->addSubMatrix(A, 0,0);
			

}

//assemblo il RHS

void ElastProblem::assembleRHS(LinearSystem* sys, std::string where)
{
	
		scalarVectorPtr_Type  source;
		source.reset(new scalarVector_Type (M_DispFEM.nb_dof()));
		gmm::clear(*source);
		
		bulkLoad(source, M_Bulk, M_DispFEM, M_CoeffFEM, M_intMethod, M_time->time());
		//givenFluidP(source, M_Bulk, M_DispFEM, M_StressScalarFEM, M_intMethod);
			
		M_Sys->addSubVector(source,0);

		
		scalarVectorPtr_Type  BCvec;
		BCvec.reset(new scalarVector_Type (M_DispFEM.nb_dof()));
		gmm::clear(*BCvec);
		//essentialWNitscheRHSVec( BCvec, M_Bulk, M_time->time(), &M_BC, M_DispFEM, M_CoeffFEM, M_intMethod);
		M_Sys->addSubVector(BCvec, 0);		
	
		gmm::clear(*BCvec);
	
		stressRHS( BCvec, M_Bulk, M_time->time(), &M_BC, M_DispFEM, M_CoeffFEM, M_intMethod);
		M_Sys->addSubVector(BCvec, 0);	
	
	
}

void ElastProblem::enforceStrongBC(bool firstTime)
{
    if (firstTime)
    {
    	for ( size_type bndID = 0; bndID < M_BC.getDiriBD().size(); bndID++ )
   	{
		dal::bit_vector quali=M_DispFEM.getFEM()->dof_on_region( M_BC.getDiriBD()[bndID]);
		
		for(dal::bv_visitor i(quali); !i.finished(); ++i)
		{
			M_rowsStrongBC.push_back(i);
			M_rowsStrongBCFlags.push_back(bndID);
		}
	}

		for (int i=0;i<M_rowsStrongBC.size(); ++++++i)
		{

			size_type ii=M_rowsStrongBC[i];
			
			bgeot::base_node where=M_DispFEM.getFEM()->point_of_basic_dof(ii);

			M_Sys->setNullRow(ii);
			M_Sys->setMatrixValue(ii,ii,1);
			scalar_type value= (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0]+
				M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0]*M_time->time());
			M_Sys->setRHSValue(ii,value);
			ii=M_rowsStrongBC[i+1];
			M_Sys->setNullRow(ii);
			M_Sys->setMatrixValue(ii,ii,1);
			value= (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1]+
				M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1]*M_time->time());
			M_Sys->setRHSValue(ii,value);
    
			ii=M_rowsStrongBC[i+2];
			M_Sys->setNullRow(ii);
			M_Sys->setMatrixValue(ii,ii,1);
			value= (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2]+
				M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2]*M_time->time());
			M_Sys->setRHSValue(ii,value);
		}
    	
    }
    else
   {
    	for (int i=0;i<M_rowsStrongBC.size(); ++++++i)
		{
			size_type ii=M_rowsStrongBC[i];
			bgeot::base_node where=M_DispFEM.getFEM()->point_of_basic_dof(ii);
			
			scalar_type value= (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0]+
				M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[0]*M_time->time());
			M_Sys->setRHSValue(ii,value);

			ii=M_rowsStrongBC[i+1];
			value= (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1]+
				M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[1]*M_time->time());
			M_Sys->setRHSValue(ii,value);
			ii=M_rowsStrongBC[i+2];
			value= (M_BC.BCDiriVec(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2]+
				M_BC.BCDiriVel(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]], M_time->time())[2]*M_time->time());
			M_Sys->setRHSValue(ii,value);
		}
    }

}

// risolve il sistema e estrae la soluzione

void ElastProblem::solve()
{
       	M_Sys->solve();

       	M_DispSol.reset(new scalarVector_Type (M_DispFEM.nb_dof()));
	gmm::clear(*M_DispSol);
	
       	M_Sys->extractSubVector(M_DispSol, 0, "sol");
}



void ElastProblem::extractSol(scalarVectorPtr_Type sol)
{
	M_DispSol.reset(new scalarVector_Type (M_DispFEM.nb_dof()));

	gmm::clear(*M_DispSol);
	gmm::copy ( gmm::sub_vector (*sol, gmm::sub_interval (0,  M_DispFEM.nb_dof()) ), *M_DispSol);
}

void ElastProblem::updateSol()
{
	gmm::clear(*M_DispSolOld);
	gmm::copy(*M_DispSol, *M_DispSolOld);

}

scalarVector_Type ElastProblem::getNormalStressOnFault()
{
	scalarVector_Type provv(M_StressScalarFEM.nb_dof("base"));
	for (size_type i=0; i<M_StressScalarFEM.nb_dof("extended"); ++i)
	{
		provv[M_StressScalarFEM.getExt(i)]=0.5*((*M_NormalStressScalarSol)[i+M_StressScalarFEM.nb_dof("base")]+(*M_NormalStressScalarSol)[M_StressScalarFEM.getExt(i)]);
	}
	return provv;
}


//
scalar_type ElastProblem::computeError(std::string what, scalar_type time)
{
    scalar_type error = 0.0;
    
    if (what == "Displacement" || what == "all")
    {
        // Create local error vector
        scalarVector_Type loc_err(M_DispFEM.nb_dof());
        
        // Compute pointwise error
        size_type qdim = M_DispFEM.get_qdim(); // Get dimension of vector field
        for (size_type i = 0; i < M_DispFEM.nb_dof(); ++i)
        {
            // For vector fields, dofs are interleaved (x0,y0,z0,x1,y1,z1,...)
            size_type comp = i % qdim; // Component (0=x, 1=y, 2=z)
            size_type base_dof = i / qdim; // Base dof index
            
            bgeot::base_node where = M_DispFEM.getFEM()->point_of_basic_dof(base_dof);
            base_small_vector exact_val = M_Bulk->getElastData()->uEx(where, time);
            
            loc_err[i] = (*M_DispSol)[i] - exact_val[comp];
        }
        
        // Compute L2 norm of error
        error = L2Norm_Elast(M_dispMassMatrix, loc_err, M_Bulk, *(M_DispFEM.getFEM()), M_intMethod, -1);
        std::cout << "at time " << time << " error displacement " << error << std::endl;
    }
    
    return error;
}


void ElastProblem::exportVtk(std::string folder, std::string what, int frame)
{
	std::cout << "export elast solution"<<std::endl;
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
	
	 	
	if (what=="Displacement" || what=="all")
	{

		getfem::vtk_export exp(folder + "Displacement_bulk"+ frameNum + ".vtk" );
  		

			exp.exporting( *(M_DispFEM.getFEM()));
			std::vector<scalar_type> dispCut(M_DispFEM.getFEM()->nb_dof(),0.0);
			gmm::copy(*M_DispSol, dispCut);
			exp.write_mesh();

			exp.write_point_data( *(M_DispFEM.getFEM()), dispCut, "u");

		getfem::vtk_export expE(folder + "Displacement_exact"+ frameNum + ".vtk" );

			expE.exporting( *(M_DispFEM.getFEM()));
			gmm::clear(dispCut);
			for (int i=0; i<dispCut.size();++++++i )
			{	
				bgeot::base_node where=M_DispFEM.getFEM()->point_of_basic_dof(i);
				dispCut[i]= M_Bulk->getElastData()->uEx(where,M_time->time())[0] ;
				dispCut[i+1]= M_Bulk->getElastData()->uEx(where,M_time->time())[1] ;
				dispCut[i+2]= M_Bulk->getElastData()->uEx(where,M_time->time())[2] ;
			}			
			
			expE.write_mesh();

			expE.write_point_data( *(M_DispFEM.getFEM()), dispCut, "u");


		
		
	}


}
