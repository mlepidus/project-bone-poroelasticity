// Core_Modern.h - Proposed modernization
#ifndef CORE_MODERN_H
#define CORE_MODERN_H

#include "Core.h"

#include <memory>
#include <vector>
#include <string>
#include <optional>
#include <variant>
#include <chrono>

// Use std::shared_ptr instead of boost::shared_ptr
template<typename T>
using SharedPtr = std::shared_ptr<T>;

template<typename T>
using UniquePtr = std::unique_ptr<T>;

// Type aliases for clarity
using ScalarType = double;
using SizeType = std::size_t;
using IndexType = std::ptrdiff_t;

// Sparse matrix types
using SparseVector = gmm::rsvector<ScalarType>;
using SparseMatrix = gmm::row_matrix<SparseVector>;
using SparseMatrixPtr = SharedPtr<SparseMatrix>;

using ScalarVector = std::vector<ScalarType>;
using ScalarVectorPtr = SharedPtr<ScalarVector>;

// Result type for operations that can fail
template<typename T>
using Result = std::variant<T, std::string>;  // Value or error message

#endif // CORE_MODERN_H