#ifndef BART_SRC_DATA_MATRICES_H_
#define BART_SRC_DATA_MATRICES_H_

#include <memory>
#include <unordered_map>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include "matrix_parameters.h"

namespace bart {

namespace data {

typedef dealii::PETScWrappers::MPI::SparseMatrix MPISparseMatrix;
typedef dealii::ConstraintMatrix                 ConstraintMatrix;

struct Matrices {
  std::unordered_map<int, std::shared_ptr<MPISparseMatrix>> system_lhs;
  std::unordered_map<int, std::shared_ptr<ConstraintMatrix>> constraints;
};

std::shared_ptr<MPISparseMatrix> BuildMatrix(MatrixParameters &parameters);

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_MATRICES_H_
