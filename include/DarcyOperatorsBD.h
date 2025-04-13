#ifndef DARCYOPERATORSBD_H
#define DARCYOPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"


void essentialWNitsche( sparseMatrixPtr_Type M,
               Bulk* medium, BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im);   //condizione al bordo di Dirichlet imposta con penalizzazione - pezzo LHS
void essentialWNitscheRHS( scalarVectorPtr_Type V,
               Bulk* medium, BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im);   //condizione al bordo di Dirichlet imposta con penalizzazione - pezzo RHS
void naturalRHS( scalarVectorPtr_Type V,
               Bulk* medium, BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im, scalar_type time);   //condizione al bordo naturale

#endif
