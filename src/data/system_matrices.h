#ifndef BART_SRC_DATA_MATRICES_H_
#define BART_SRC_DATA_MATRICES_H_

#include <memory>
#include <unordered_map>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

namespace bart {

namespace data {

struct Matrices {
  typedef dealii::PETScWrappers::MPI::SparseMatrix MPISparseMatrix;
  typedef dealii::ConstraintMatrix                 ConstraintMatrix;

  std::unordered_map<int, std::unique_ptr<MPISparseMatrix>> system_lhs;
  std::unordered_map<int, std::unique_ptr<ConstraintMatrix>> constraints;
};

} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_MATRICES_H_
