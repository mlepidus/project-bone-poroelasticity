// ============================================================================
// Elasticity Operator Headers
// ============================================================================

// ElastOperatorsBulk.h - Assembly operators for bulk elasticity
#ifndef ELASTOPERATORSBULK_H
#define ELASTOPERATORSBULK_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

/**
 * @brief Assemble elastic stiffness matrix
 * @param M Output stiffness matrix
 * @param medium Bulk domain pointer
 * @param femV Displacement finite element space
 * @param femP Coefficient finite element space
 * @param im Integration method
 */
void stiffElast(sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP,
                getfem::mesh_im& im);

/**
 * @brief Assemble volumetric body force RHS
 * @param V Output RHS vector
 * @param medium Bulk domain pointer
 * @param FemD Displacement finite element space
 * @param FemC Coefficient finite element space
 * @param im Integration method
 * @param t Current time
 */
void bulkLoad(scalarVectorPtr_Type V, Bulk* medium, FEM& FemD, FEM& FemC,
              getfem::mesh_im& im, scalar_type t);

/**
 * @brief Apply prescribed fluid pressure field (decoupled/split problem)
 * @param V Output RHS vector
 * @param medium Bulk domain pointer
 * @param FemD Displacement finite element space
 * @param FemC Coefficient finite element space
 * @param im Integration method
 */
void givenFluidP(scalarVectorPtr_Type V, Bulk* medium, FEM& FemD, FEM& FemC,
                 getfem::mesh_im& im);

/// Apply fluid pressure from solution vector
void givenFluidP(scalarVectorPtr_Type V, scalarVectorPtr_Type pressure,
                 Bulk* medium, FEM& FemD, FEM& FemC, getfem::mesh_im& im);

/**
 * @brief Assemble poroelastic coupling matrix (fully coupled problem)
 * @param M Output coupling matrix (pressure affects stress)
 * @param FemD Displacement finite element space
 * @param FemC Coefficient finite element space
 * @param im Integration method
 */
void matrixFluidP(sparseMatrixPtr_Type M, FEM& FemD, FEM& FemC,
                  getfem::mesh_im& im);

/**
 * @brief Assemble mass matrix for displacement field
 * @param M Output mass matrix
 * @param FemD Displacement finite element space
 * @param im Integration method
 */
void massMatrix(sparseMatrixPtr_Type M,  FEM& FemD,
                getfem::mesh_im& im);

/**
 * @brief Compute L2 norm of elastic displacement error
 * @param M Mass matrix for integration
 * @param V Solution vector
 * @return L2 norm value
 */
scalar_type L2Norm_Elast(sparseMatrixPtr_Type M, scalarVector_Type V);

#endif // ELASTOPERATORSBULK_H
