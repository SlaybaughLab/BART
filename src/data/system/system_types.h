#ifndef BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_
#define BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>

namespace bart {

namespace data {

namespace system {

//! Group number for rhs and lhs
using GroupNumber = int;
//! Angle index for rhs and lhs
using AngleIndex = int;
//! Index used to store and access rhs vectors and lhs matrices
using Index = std::pair<GroupNumber, AngleIndex>;

//! Sparse MPI matrix used for left-hand-side matrices
using MPISparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
//! Sparse MPI vector used for ride-hand-side vectors
using MPIVector = dealii::PETScWrappers::MPI::Vector;

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_