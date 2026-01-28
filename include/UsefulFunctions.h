// ============================================================================
// UsefulFunctions.h - General utility functions for mesh and solution export
// ============================================================================
#ifndef USEFULFUNCTIONS_H
#define USEFULFUNCTIONS_H

#include "Core.h"
#include "FEM.h"

// Export solution field to VTK format for visualization
void exportSolution(const std::string& fileName,
                    const std::string& solutionName,
                    const getfem::mesh_fem& meshFEM,
                    const scalarVector_Type& solution);

// Export solution as cell-centered data (discontinuous representation)
void exportSolutionInCell(const std::string& fileName,
                          const std::string& solutionName,
                          const getfem::mesh_fem& meshFEM,
                          const scalarVector_Type& solution);

// Export mesh structure to VTK format
void exportMesh(const std::string& fileName, const getfem::mesh& mesh);

// Convert mass matrix to lumped (diagonal) form for explicit time integration
void massLumping(sparseMatrix_Type& matrix);

// Convert GetFem++ bit_vector to standard vector for easier manipulation
void fromBitVectorToStdVector(dal::bit_vector& bitVector,
                              std::vector<size_type>& stdVector);

// Utility function to convert integer to character representation
char intToChar(const size_type& integer);

// Generate region sign string from level set values
std::string regionSigns(const scalarVector_Type& levelSetValue);

// Compute Euclidean distance between two 2D points
scalar_type pointDistance(const scalar_type& x0, const scalar_type& x1,
                         const scalar_type& y0, const scalar_type& y1);

// Determine boolean operation for sub-region definition with level sets
std::string getOperation(const std::string& subRegion,
                        const sizeVector_Type& levelSets);

// Compare sign patterns for region matching
std::pair<std::string, size_type> comparaSegni(const std::string& region,
                                               const scalarVector_Type& signs);

// Check if a point lies inside a triangular element
bool isInTriangle(const getfem::mesh& mesh,
                 const size_type& elementID,
                 const base_node& node,
                 const scalar_type& toll = 1e-7);

// Check if a point lies inside a triangle defined by three vertices
bool isInTriangle(const base_node& v1, const base_node& v2, const base_node& v3,
                 const base_node& node,
                 const scalar_type& toll = 1e-7);

// Compute intersection point between two line segments
bool intersectSegments(const base_node& a1, const base_node& a2,
                      const base_node& b1, const base_node& b2,
                      base_node& sol);

// Compute intersection between a line segment and a triangle
bool intersectSegmentTriangle(const base_node& a1, const base_node& a2,
                             const base_node& v1, const base_node& v2, const base_node& v3,
                             base_node& s1, base_node& s2);

#endif // USEFULFUNCTIONS_H
