#ifndef BART_SRC_DATA_MATRIX_PARAMETERS_
#define BART_SRC_DATA_MATRIX_PARAMETERS_

#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

/*! This struct holds the parameters necessary for creation of a parallel
 * PETSc matrix. */

namespace bart {

namespace data {

typedef dealii::PETScWrappers::MPI::SparseMatrix MPISparseMatrix;

struct MatrixParameters {
  dealii::IndexSet rows;
  dealii::IndexSet columns;
  dealii::DynamicSparsityPattern sparsity_pattern;
};

std::shared_ptr<MPISparseMatrix> BuildMatrix(MatrixParameters &parameters); 

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_MATRIX_PARAMETERS_
