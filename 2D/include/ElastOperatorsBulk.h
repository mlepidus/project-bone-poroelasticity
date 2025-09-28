#ifndef ELASTOPERATORSBULK_H
#define ELASTOPERATORSBULK_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

//operatori elasticità nel bulk

void stiffElast ( sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP, getfem::mesh_im& im);  //matrice di stiffness
void bulkLoad (scalarVectorPtr_Type V, Bulk* medium,  FEM& FemD, FEM& FemC, getfem::mesh_im& im, scalar_type t);  //termine forzante volumetrico 
void givenFluidP(scalarVectorPtr_Type V, Bulk* medium,  FEM& FemD, FEM& FemC, getfem::mesh_im& im);  //pressione del fluido imposta (per problema disaccoppiato o splittato)
void givenFluidP(scalarVectorPtr_Type V, scalarVectorPtr_Type pressure, Bulk* medium,  FEM& FemD, FEM& FemC, getfem::mesh_im& im);  //pressione del fluido imposta (per problema disaccoppiato o splittato)
void matrixFluidP(sparseMatrixPtr_Type M,  Bulk* medium,  FEM& FemD, FEM& FemC, getfem::mesh_im& im);  //matrice di accoppiamento poroelastico per problema fully coupled  
void massMatrix(sparseMatrixPtr_Type M,  Bulk* medium,  FEM& FemD, getfem::mesh_im& im);  //matrice di massa per i fem del displacement

#endif
