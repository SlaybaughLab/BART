#include "solver/eigenvalue/krylov_schur_eigenvalue_solver.hpp"

#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/petsc_vector.h>

namespace bart::solver::eigenvalue {

auto KrylovSchurEigenvalueSolver::SpectralRadius(const SpectralRadiusI::MatrixBase& base)
-> std::pair<double, std::vector<double>> {
  dealii::SolverControl solver_control(1000, 1e-6);
  dealii::SLEPcWrappers::SolverKrylovSchur solver(solver_control);

  std::vector<double> eigenvalue(1);
  std::vector<dealii::PETScWrappers::MPI::Vector> eigenvector(1);
  eigenvector.at(0).reinit(MPI_COMM_WORLD, base.m(), base.m());

  solver.solve(base, eigenvalue, eigenvector);

  std::vector<double> return_eigenvector(eigenvector.at(0).size());
  for (int i = 0; i < return_eigenvector.size(); ++i)
    return_eigenvector.at(i) = eigenvector.at(0)[i];

  return {eigenvalue.at(0), return_eigenvector};
}

} // namespace bart::solver::eigenvalue
