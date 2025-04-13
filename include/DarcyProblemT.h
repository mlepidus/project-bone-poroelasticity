#ifndef DARCYPROBLEMT_H
#define DARCYPROBLEMT_H

#include "LinearSystem.h"
#include "DarcyOperators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "TimeLoop.h"

//classe per il problema di darcy tempo dipendente. rispetto al problema stazionario contiene anche un puntatore al ciclo in tempo

class DarcyProblemT
{
public: 
	 DarcyProblemT ( const GetPot& dataFile, Bulk* bulk=NULL);

	 inline void setTimePointer(TimeLoop* timePtr)
	 {	
		M_timeLoop=timePtr;
		M_dt=timePtr->dt();
	 }
	 
	 FEM* getFEM(const std::string where="bulk", const std::string variable="Pressure");
	 
	 void initialize();
	 
         void addToSys(LinearSystem* sys);
         
         void assembleMatrix(LinearSystem* sys, std::string where);
         
         void assembleRHS(LinearSystem* sys, std::string where);
         
         void solve();
         
         void updateSol();
                
         void exportVtk(std::string folder="./vtk", std::string what="all", int frame=-1);

	 scalar_type computeError(std::string what="all", scalar_type time=0);
         
         size_type getNDOF(std::string variable="all");
         
         void extractSol(scalarVectorPtr_Type sol);
          
        
 private:

    scalar_type M_dt;
    TimeLoop* M_timeLoop;
    Bulk* M_Bulk;
    BC M_BC;
    FEM M_PressureFEM;

    FEM M_CoeffFEM;

    FEM M_VelocityFEM;

    FEM M_visualizationVFEM;
    
    LinearSystem* M_Sys;
    
    scalarVectorPtr_Type M_pressureSol;
    scalarVectorPtr_Type M_pressureSolOld;
    scalarVectorPtr_Type M_pressureSolIni;
    scalarVectorPtr_Type M_velocitySol;

    sparseMatrixPtr_Type A22;


    size_type M_nbTotDOF;
    size_type M_nbTotBulkDOF;
    
    getfem::mesh_im M_intMethod;
    


    mutable LifeV::Parser M_parser;
};


#endif

