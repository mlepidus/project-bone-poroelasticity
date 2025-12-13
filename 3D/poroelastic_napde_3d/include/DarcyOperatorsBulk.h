// ============================================================================
// Darcy Operator Headers
// ============================================================================

// DarcyOperatorsBulk.h - Assembly operators for bulk Darcy flow
#ifndef DARCYOPERATORSBULK_H
#define DARCYOPERATORSBULK_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

/**
 * @brief Assemble H(div) mass matrix for velocity field
 * @param M Output mass matrix A in mixed formulation
 * @param medium Bulk domain pointer
 * @param femV Velocity finite element space
 * @param femP Pressure finite element space
 * @param im Integration method
 */
void massHdiv(sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP,
              getfem::mesh_im& im);

/**
 * @brief Assemble L2 mass matrix for pressure
 * @param M Output mass matrix (scaled by time step)
 * @param medium Bulk domain pointer
 * @param femV Velocity finite element space
 * @param femP Pressure finite element space
 * @param im Integration method
 * @param dt Time step size
 */
void massL2(sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP,
            getfem::mesh_im& im, scalar_type dt);

/// Local mass matrix assembly for control volume icv
void massL2Local(sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP,
                 getfem::mesh_im& im, scalar_type dt, size_type icv);

/// Mass matrix assembly on specific mesh region
void massL2(sparseMatrixPtr_Type M, Bulk* medium, getfem::mesh_fem& femV,
            getfem::mesh_fem& femP, getfem::mesh_im& im, scalar_type dt, int region);

/// Compute L2 norm of scalar field
scalar_type L2Norm(scalarVector_Type V, Bulk* medium, FEM& femP, getfem::mesh_im& im);

/// Compute L2 norm of vector field
scalar_type L2Norm(scalarVector_Type V, Bulk* medium, FEM& femV, FEM& femC,
                   getfem::mesh_im& im);

/// Compute L2 norm using pointer to vector
scalar_type L2Norm(scalarVectorPtr_Type V, Bulk* medium, FEM& femP,
                   getfem::mesh_im& im);

/// Compute weighted L2 norm with mass matrix M
scalar_type L2Norm(sparseMatrixPtr_Type M, scalarVector_Type V, Bulk* medium,
                   getfem::mesh_fem& femP, getfem::mesh_im& im, int region = -1);

/**
 * @brief Assemble divergence operator B in mixed formulation
 * @param M Output divergence matrix (couples velocity to pressure)
 * @param medium Bulk domain pointer
 * @param femV Velocity finite element space
 * @param femP Pressure finite element space
 * @param im Integration method
 */
void divHdiv(sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP,
             getfem::mesh_im& im);

/// Assemble scalar source term for divergence equation RHS
void scalarSource(scalarVectorPtr_Type V, Bulk* medium, FEM& femP, FEM& femC,
                  getfem::mesh_im& im, const scalar_type t);

/// Assemble vector source term (typically gravity)
void vectorSource(scalarVectorPtr_Type V, Bulk* medium, FEM& femV, FEM& femC,
                  getfem::mesh_im& im);

/// Assemble vector source term from data field
void vectorSource(scalarVectorPtr_Type V, Bulk* medium, FEM& femV, FEM& femC,
                  scalarVectorPtr_Type data, getfem::mesh_im& im);

#endif // DARCYOPERATORSBULK_H