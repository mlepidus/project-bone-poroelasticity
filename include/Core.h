// ============================================================================
// Core.h - Foundation header with all essential includes and type definitions
// ============================================================================
#ifndef CORE_H
#define CORE_H

// Standard library includes
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>

// Boost library
#include <boost/shared_ptr.hpp>

// GetFem++ core includes
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_config.h>
#include <getfem/getfem_assembling.h>
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>

// GetFem++ mesh structures
#include <getfem/bgeot_mesh.h>

// GetFem++ level set and XFEM support
#include <getfem/getfem_mesh_im_level_set.h>
#include <getfem/getfem_mesh_fem_level_set.h>
#include <getfem/getfem_mesh_fem_product.h>
#include <getfem/getfem_mesh_fem_global_function.h>

// GMM++ linear algebra library
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_MUMPS_interface.h>
#include <gmm/gmm_superlu_interface.h>

// GetPot for parameter file parsing
#include "GetPot"

// Import common GetFem++ types into current namespace
using bgeot::base_small_vector;  // Specialized vector class for small dimensions (dim < 16)
using bgeot::base_node;          // Geometric node representation (derived from base_small_vector)
using bgeot::scalar_type;        // Floating point type (typically double)
using bgeot::size_type;          // Unsigned integer type for sizes and indices

// ============================================================================
// Type Definitions - Sparse matrix and vector types using GMM++
// ============================================================================

// Sparse vector type using row-sparse storage
typedef gmm::rsvector<scalar_type> sparseVector_Type;

// Sparse matrix types
typedef gmm::row_matrix<sparseVector_Type> sparseMatrix_Type;
typedef boost::shared_ptr<sparseMatrix_Type> sparseMatrixPtr_Type;
typedef std::vector<sparseMatrix_Type> sparseMatrixContainer_Type;
typedef std::vector<sparseMatrixPtr_Type> sparseMatrixPtrContainer_Type;

// Dense vector types
typedef std::vector<scalar_type> scalarVector_Type;
typedef boost::shared_ptr<scalarVector_Type> scalarVectorPtr_Type;
typedef std::vector<scalarVector_Type> scalarVectorContainer_Type;
typedef boost::shared_ptr<scalarVectorPtr_Type> scalarVectorContainerPtr_Type;
typedef std::vector<scalarVectorPtr_Type> scalarVectorPtrContainer_Type;

// Size/index vector types
typedef std::vector<size_type> sizeVector_Type;
typedef boost::shared_ptr<sizeVector_Type> sizeVectorPtr_Type;
typedef std::vector<sizeVector_Type> sizeVectorContainer_Type;
typedef std::vector<sizeVectorPtr_Type> sizeVectorPtrContainer_Type;

// Pair and string container types
typedef std::pair<size_type, size_type> pairSize_Type;
typedef std::vector<pairSize_Type> pairSizeVector_Type;
typedef std::vector<pairSizeVector_Type> pairSizeVectorContainer_Type;
typedef std::vector<std::string> stringContainer_Type;

// GetFem++ mesh and finite element types
typedef getfem::mesh_fem GFMeshFEM_Type;
typedef boost::shared_ptr<GFMeshFEM_Type> GFMeshFEMPtr_Type;
typedef std::vector<GFMeshFEMPtr_Type> GFMeshFEMPtrContainer_Type;

// GetFem++ level set types
typedef getfem::mesh_level_set GFMeshLevelSet_Type;
typedef boost::shared_ptr<GFMeshLevelSet_Type> GFMeshLevelSetPtr_Type;
typedef std::vector<GFMeshLevelSetPtr_Type> GFMeshLevelSetPtrContainer_Type;

typedef getfem::level_set GFLevelSet_Type;
typedef std::vector<GFLevelSet_Type> GFLevelSetContainer_Type;
typedef boost::shared_ptr<GFLevelSet_Type> GFLevelSetPtr_Type;
typedef std::vector<GFLevelSetPtr_Type> GFLevelSetPtrContainer_Type;

typedef getfem::mesh_im_level_set GFIntegrationMethodLevelSet_Type;
typedef boost::shared_ptr<GFIntegrationMethodLevelSet_Type> GFIntegrationMethodLevelSetPtr_Type;

// LifeV compatibility layer - type aliases for legacy code
namespace LifeV {
    typedef scalar_type Double;
    typedef size_type UInt;
    typedef size_type ID;
    typedef int Int;
}

#endif // CORE_H