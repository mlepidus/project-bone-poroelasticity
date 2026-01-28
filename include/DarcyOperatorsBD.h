// DarcyOperatorsBD.h - Boundary operators for Darcy flow
#ifndef DARCYOPERATORSBD_H
#define DARCYOPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

/**
 * @brief Assemble Nitsche penalty term for essential BC (LHS contribution)
 * @param M Output matrix (penalty terms)
 * @param medium Bulk domain pointer
 * @param bcPtr Boundary condition handler
 * @param femV Velocity finite element space
 * @param femP Pressure finite element space
 * @param im Integration method
 */
void essentialWNitsche(sparseMatrixPtr_Type M, Bulk* medium, BC* bcPtr,
                       FEM& femV, FEM& femP, getfem::mesh_im& im);

/**
 * @brief Assemble Nitsche penalty term for essential BC (RHS contribution)
 * @param V Output RHS vector
 * @param medium Bulk domain pointer
 * @param bcPtr Boundary condition handler
 * @param femV Velocity finite element space
 * @param femP Pressure finite element space
 * @param im Integration method
 */
void essentialWNitscheRHS(scalarVectorPtr_Type V, Bulk* medium, BC* bcPtr,
                          FEM& femV, FEM& femP, getfem::mesh_im& im);

/**
 * @brief Assemble natural (Neumann) boundary condition RHS
 * @param V Output RHS vector
 * @param medium Bulk domain pointer
 * @param bcPtr Boundary condition handler
 * @param femV Velocity finite element space
 * @param femP Pressure finite element space
 * @param im Integration method
 * @param time Current time for time-dependent BC
 */
void naturalRHS(scalarVectorPtr_Type V, Bulk* medium, BC* bcPtr,
                FEM& femV, FEM& femP, getfem::mesh_im& im, scalar_type time);

#endif // DARCYOPERATORSBD_H