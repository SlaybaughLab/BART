#ifndef BART_SRC_SOLVER_GMRES_H_
#define BART_SRC_SOLVER_GMRES_H_

#include <memory>

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_vector_base.h>

#include "solver/linear_i.h"

namespace bart {

namespace solver {

class GMRES : public LinearI {
 public:
  GMRES(int max_iterations = 100, double convergence_tolerance = 1e-10);
  ~GMRES() = default;

  void Solve(dealii::PETScWrappers::MatrixBase *A,
             dealii::PETScWrappers::VectorBase *x,
             dealii::PETScWrappers::VectorBase *b,
             dealii::PETScWrappers::PreconditionerBase *preconditioner) override;
  int max_iterations() const { return solver_control_.max_steps(); };
  double convergence_tolerance() const { return solver_control_.tolerance(); };

  const dealii::SolverControl& solver_control() const { return solver_control_;};

 private:
  dealii::SolverControl solver_control_;
};

} // namespace solver

} // namespace bart

#endif // BART_SRC_SOLVER_GMRES_H_