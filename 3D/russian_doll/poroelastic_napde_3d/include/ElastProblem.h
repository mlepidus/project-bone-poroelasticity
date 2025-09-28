//classe per la gestione del problema dell'elasticità

#ifndef ELASTPROBLEM_H
#define ELASTPROBLEM_H

#include "LinearSystem.h"
#include "ElastOperators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "TimeLoop.h"


class ElastProblem
{
public: 
	 ElastProblem ( const GetPot& dataFile, Bulk* bulk=NULL);

	inline void setTimePointer(TimeLoop* timePtr)  //anche se quasi statico le BC, i carichi non lo sono, quindi devo conoscere il tempo
	{
		M_time=timePtr;
	}
	 
	 FEM* getFEM(const std::string where="bulk", const std::string variable="Pressure");




	 inline BC* getBC()
	 {
		return &M_BC;
      }
	 
         void addToSys(LinearSystem* sys); 
         
         void assembleMatrix(LinearSystem* sys, std::string where);
         
         void assembleRHS(LinearSystem* sys, std::string where);

	void enforceStrongBC(bool firstTime);

   scalar_type computeError(std::string what, scalar_type time);
	 
   inline LinearSystem* getSys()
	 {
		return M_Sys;
	 }
         
         void solve();
         
         void updateSol();
         
         void extractSol(scalarVectorPtr_Type sol);
                
         void exportVtk(std::string folder="./vtk", std::string what="all", int frame=-1);
         
         size_type getNDOF(std::string variable="all");
         
         scalarVectorPtr_Type getOldSol()
         {
         	return M_DispSolOld;
         }
         
         void initialize();
         
         scalarVector_Type getNormalStressOnFault();



	 inline scalar_type frictionCoeff()
	 {
		return M_staticMu;
	 }
         
 private:

    TimeLoop* M_time;
    Bulk* M_Bulk;

    BC M_BC;
    FEM M_DispFEM;
    FEM  M_DispScalarFEM;
    FEM M_CoeffFEM;
    FEM M_StressFEM;
    FEM M_StressScalarFEM;

    LinearSystem* M_Sys;

    scalarVectorPtr_Type M_DispSol;
    scalarVectorPtr_Type M_DispSolOld;
    scalarVectorPtr_Type M_slipVel;
    scalarVectorPtr_Type M_NormalStressSol;
    scalarVectorPtr_Type M_NormalStressScalarSol;
    scalarVectorPtr_Type M_TangentStressScalarSol;
    scalarVectorPtr_Type M_ST;
    scalarVectorPtr_Type M_Slip;
    scalarVectorPtr_Type M_SlipOld;
    sparseMatrixPtr_Type M_normalStressMat;
    sparseMatrixPtr_Type M_normalStressScalarMat;
    sparseMatrixPtr_Type M_tangentStressScalarMat;
    sparseMatrixPtr_Type M_massMatrixVector;
    sparseMatrixPtr_Type M_massMatrix;
    sparseMatrixPtr_Type M_massP2onFault;
    sparseMatrixPtr_Type M_maskUzawa;
    sparseMatrixPtr_Type M_glueMatrix;
    sparseMatrixPtr_Type M_normalGlueMatrix;
    sparseMatrixPtr_Type M_P02P2interp;

    sparseMatrixPtr_Type M_dispMassMatrix; //

    size_type M_nbTotDOF;
    size_type M_nbTotBulkDOF;
    std::vector<size_type> M_rowsStrongBC;
    std::vector<size_type> M_rowsStrongBCFlags;

    scalar_type M_staticMu;
   
    bool M_precomputeInterpolation;
    
    getfem::mesh_im M_intMethod;

    bool M_IsNitSym;
    bool M_IsNitCons;

   // mutable LifeV::Parser M_parser;
};


#endif

