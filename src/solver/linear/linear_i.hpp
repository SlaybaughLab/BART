#ifndef BART_SOLVER_LINEAR_I_H_
#define BART_SOLVER_LINEAR_I_H_

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_vector_base.h>

namespace bart::solver::linear {
/*! \brief Linear solver class.
 *
 * This class is all solvers that solve an equation in the form \f$Ax = b\f$.
 *
 */
class LinearI {
 public:
  virtual ~LinearI() = default;
  virtual void Solve(
      dealii::PETScWrappers::MatrixBase *A,
      dealii::PETScWrappers::VectorBase *x,
      dealii::PETScWrappers::VectorBase *b,
      dealii::PETScWrappers::PreconditionerBase *preconditioner) = 0;
};

} // namespace bart::solver::linear

#endif // BART_SOLVER_LINEAR_I_H_