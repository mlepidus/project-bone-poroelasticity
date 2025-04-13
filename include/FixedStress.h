#ifndef FIXEDSTRESS_H
#define FIXEDSTRESS_H

#include "LinearSystem.h"
//#include "CouplingEDOperators.h"
#include "ElastProblem.h"
#include "ElastOperators.h"
#include "DarcyProblemT.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "TimeLoop.h"

//classe per la gestione del loop iterativo fixed stress per il problema fluido+meccanico

class FixedStress
{
public: 
	 FixedStress ( const GetPot& dataFile, Bulk* bulk=NULL,  TimeLoop* time=NULL );
	 
	 void inline addElastPB(ElastProblem* elast)   //prendo il puntatore a pb elasticità
	 {
	 	M_ElastPB=elast;
		M_ElastPB->setTimePointer(M_time);
	 }
	 
	 void inline addDarcyPB(DarcyProblemT* darcy) //prendo il puntatore a pb fluido
	 {
	 	M_DarcyPB=darcy;	
		M_DarcyPB->setTimePointer(M_time);
	 }

         void initialize();
	
	 void assembleMatrices();   //chiamo gli assembler dei due

         void solve();

	 void updateSol();

	 void solveFluidStep();

	 void solveMechStep();

	 scalar_type stoppingCriterion();
                
         void exportVtk(std::string folder="./vtk", std::string what="all", int frame=-1);
         
 private:

   
    Bulk* M_Bulk;
    
    TimeLoop* M_time;
    
    LinearSystem* M_elastSys;
    LinearSystem* M_darcySys;
    
    ElastProblem* M_ElastPB;
    DarcyProblemT* M_DarcyPB;

    scalarVectorPtr_Type M_solEl;  //n+1,i
    scalarVectorPtr_Type M_solDa;

    scalarVectorPtr_Type M_solElPre;  //n+1,i-1
    scalarVectorPtr_Type M_solDaPre;

    scalarVectorPtr_Type M_solElOld;  //n 
    scalarVectorPtr_Type M_solDaOld;
    
    sparseMatrixPtr_Type M_PressureStress;
    sparseMatrixPtr_Type M_3D2D;
    
    getfem::mesh_im M_intMethod;
    
    size_type M_nbTotDOF;
 
    std::string M_betaString;
    std::string M_contactString;
    std::vector<scalar_type> M_beta;
    scalar_type M_toll;
    scalar_type M_ContactToll;
    scalar_type M_ContactMaxIt;
    size_type M_maxit;

    sparseMatrixPtr_Type M_massMatrix;
    sparseMatrixPtr_Type M_massMatrixF;

    bool M_isCoupled;

    mutable LifeV::Parser M_parser;
};


#endif

