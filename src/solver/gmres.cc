#include "solver/gmres.h"

#include <deal.II/lac/petsc_solver.h>

namespace bart {

namespace solver {

GMRES::GMRES(int max_iterations, double convergence_tolerance)
    : max_iterations_(max_iterations),
      convergence_tolerance_(convergence_tolerance),
      solver_control_(max_iterations_, convergence_tolerance_){}

void GMRES::Solve(dealii::PETScWrappers::MatrixBase *A,
                  dealii::PETScWrappers::VectorBase *x,
                  dealii::PETScWrappers::VectorBase *b,
                  dealii::PETScWrappers::PreconditionerBase *preconditioner) {
  dealii::PETScWrappers::SolverGMRES solver(solver_control_, MPI_COMM_WORLD);
  solver.solve(*A, *x, *b, *preconditioner);
}

} // namespace solver

} // namespace bart