#include "solver/eigenvalue/krylov_schur_eigenvalue_solver.hpp"

#include <deal.II/base/mpi.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/petsc_vector.h>

namespace bart::solver::eigenvalue {

auto KrylovSchurEigenvalueSolver::SpectralRadius(const SpectralRadiusI::MatrixBase& base)
-> std::pair<double, std::vector<double>> {
  dealii::SolverControl solver_control(1000, 1e-6);
  dealii::SLEPcWrappers::SolverKrylovSchur solver(solver_control);

  const auto m{ base.m() };

  std::vector<double> eigenvalues(m);
  std::vector<dealii::PETScWrappers::MPI::Vector> eigenvectors(m);
  for (auto& eigenvector : eigenvectors)
    eigenvector.reinit(MPI_COMM_WORLD, m, m);

  solver.solve(base, eigenvalues, eigenvectors, m);

  auto max_element_index = std::distance(eigenvalues.begin(), std::max_element(eigenvalues.begin(), eigenvalues.end()));

  std::vector<double> return_eigenvector(m);

  for (int i = 0; i < return_eigenvector.size(); ++i)
    return_eigenvector.at(i) = eigenvectors.at(max_element_index)[i];

  return {eigenvalues.at(max_element_index), return_eigenvector};
}

} // namespace bart::solver::eigenvalue
