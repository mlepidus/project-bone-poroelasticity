#ifndef ELASTOPERATORSBD_H
#define ELASTOPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"

#include "BC.h"

//condizioni al contorno per il problema dell'elasticità

void essentialWNitscheVec( sparseMatrixPtr_Type M,
               Bulk* medium, BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im);  //penalizzazione per dirichlet - LHS
void essentialWNitscheRHSVec( scalarVectorPtr_Type V,
               Bulk* medium, scalar_type time,  BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im); //penalizzazione per dirichlet - RHS
void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium, scalar_type time, BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im);  //condizione di sforzo

#endif
