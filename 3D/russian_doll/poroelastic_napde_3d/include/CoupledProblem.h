#ifndef COUPLEDPROBLEM_H
#define COUPLEDPROBLEM_H

#include "LinearSystem.h"
#include "ElastProblem.h"
#include "DarcyProblemT.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "TimeLoop.h"

//gestione del problema strongly coupled (solo lineare ossia senza slip) - NON up to date

class CoupledProblem
{
public: 
	 CoupledProblem ( const GetPot& dataFile, Bulk* bulk=NULL,TimeLoop* time=NULL );
	 
	 void inline addElastPB(ElastProblem* elast)
	 {
	 	M_ElastPB=elast;
		M_ElastPB->setTimePointer(M_time);
	 }
	 
	 void inline addDarcyPB(DarcyProblemT* darcy)
	 {
	 	M_DarcyPB=darcy;	
		M_DarcyPB->setTimePointer(M_time);
	 }
	 
	 bgeot::base_node computeError(scalar_type t);

	 void build3D2DCoupling();
	 
         void addToSys(LinearSystem* sys);
         
         void assembleMatrix(LinearSystem* sys);
         
         void assembleRHS(LinearSystem* sys);
         
         void addSubSystems();
         void addSubSystemsRHS();
         
         void clearSubSystems();
         void clearSubSystemsRHS();
         
         void solve();
         
         void updateSol();
                
         void exportVtk(std::string folder="./vtk", std::string what="all", int frame=-1);

         // Funzione per esportare la storia
         void exportHistory(const scalarVector_Type& timestepData, const std::string& filename, size_t step);

	 void enforceStrongBC(bool firstTime);
         
 private:

   
    Bulk* M_Bulk;
    
    TimeLoop* M_time;
    
    LinearSystem* M_Sys;
    LinearSystem M_elastSys;
    LinearSystem M_darcySys;
    
    ElastProblem* M_ElastPB;
    DarcyProblemT* M_DarcyPB;
    
    scalarVectorPtr_Type M_solEl;
    scalarVectorPtr_Type M_solDa;
    
    sparseMatrixPtr_Type M_PressureStress;
    sparseMatrixPtr_Type M_3D2D;
    
    getfem::mesh_im M_intMethod;
    
    size_type M_nbTotDOF;
    //step counter
    size_type step;
 
    mutable LifeV::Parser M_parser;
};


#endif

