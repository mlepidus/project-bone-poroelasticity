#ifndef ELASTOPERATORSBD_H
#define ELASTOPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"

#include "BC.h"

// Boundary condition operators for the elasticity problem

void essentialWNitscheVec( sparseMatrixPtr_Type M,
               Bulk* medium, BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im);  //Dirichlet penality - LHS
void essentialWNitscheRHSVec( scalarVectorPtr_Type V,
               Bulk* medium, scalar_type time,  BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im); //Dirichlet penality - RHS
void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium, scalar_type time, BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im);  //stress condition - RHS

#endif
