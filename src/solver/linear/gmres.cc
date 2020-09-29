#include "gmres.h"

#include "solver/factory/solver_factory_registrar.h"

#include <deal.II/lac/petsc_solver.h>

namespace bart {

namespace solver {

factory::LinearIFactoryRegistrar<GMRES, int, double> _gmres_registration(
    LinearSolverName::kGMRES,
    [] (int max_iterations, double convergence_tolerance) {
      std::unique_ptr<LinearI> return_ptr;
      return_ptr = std::make_unique<GMRES>(max_iterations, convergence_tolerance);
      return return_ptr; }
);

GMRES::GMRES(int max_iterations, double convergence_tolerance)
    : solver_control_(max_iterations, convergence_tolerance){}

void GMRES::Solve(dealii::PETScWrappers::MatrixBase *A,
                  dealii::PETScWrappers::VectorBase *x,
                  dealii::PETScWrappers::VectorBase *b,
                  dealii::PETScWrappers::PreconditionerBase *preconditioner) {
  dealii::PETScWrappers::SolverGMRES solver(solver_control_, MPI_COMM_WORLD);
  solver.solve(*A, *x, *b, *preconditioner);
}

} // namespace solver

} // namespace bart