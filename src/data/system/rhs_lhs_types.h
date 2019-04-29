#ifndef BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_
#define BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>

namespace bart {

namespace data {

namespace system {

//! Sparse MPI matrix used for left-hand-side matrices
using MPISparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
//! Sparse MPI vector used for ride-hand-side vectors
using MPIVector = dealii::PETScWrappers::MPI::Vector;

} // namespace system

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_SYSTEM_RHS_LHS_TYPES_H_