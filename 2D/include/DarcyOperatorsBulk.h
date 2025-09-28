#ifndef DARCYOPERATORSBULK_H
#define DARCYOPERATORSBULK_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"

#include "BC.h"

// sistema lineare A B\\ B' M  (Darcy è risolto in forma mista)

void massHdiv ( sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP, getfem::mesh_im& im);  //matrice di massa per le velocità (A)
void massL2 ( sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP, getfem::mesh_im& im, scalar_type dt);  //matrice di massa, M
void massL2Local ( sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP, getfem::mesh_im& im, scalar_type dt, size_type icv);  //matrice di massa, M
void massL2 ( sparseMatrixPtr_Type M, Bulk* medium, getfem::mesh_fem& femV,  getfem::mesh_fem& femP, getfem::mesh_im& im, scalar_type dt, int region);  //matrice di massa, M
scalar_type L2Norm ( scalarVector_Type V, Bulk* medium, FEM& femP, getfem::mesh_im& im);  //matrice di massa, M
scalar_type L2Norm ( scalarVector_Type V, Bulk* medium, FEM& femV, FEM& femC,  getfem::mesh_im& im);  //matrice di massa, M
scalar_type L2Norm ( scalarVectorPtr_Type V, Bulk* medium, FEM& femP, getfem::mesh_im& im);  //matrice di massa, M

scalar_type L2Norm (sparseMatrixPtr_Type M, scalarVector_Type V, Bulk* medium, getfem::mesh_fem& femP, getfem::mesh_im& im, int region=-1);  //matrice di massa, M

void divHdiv ( sparseMatrixPtr_Type M, Bulk* medium,  FEM& femV, FEM& femP, getfem::mesh_im& im);  //matrice B
void scalarSource (scalarVectorPtr_Type V, Bulk* medium,  FEM& femP, FEM& femC, getfem::mesh_im& im, const scalar_type t);  //termine sorgente RHS per eqt della divergenza
void vectorSource (scalarVectorPtr_Type V, Bulk* medium,  FEM& femV, FEM& femC, getfem::mesh_im& im);  // termine forzante vettoriale (tipicamente gravità)
void vectorSource (scalarVectorPtr_Type V, Bulk* medium,  FEM& femV, FEM& femC, scalarVectorPtr_Type data,getfem::mesh_im& im);  // termine forzante vettoriale (tipicamente gravità)
#endif
