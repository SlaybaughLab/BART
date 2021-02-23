#include "calculator/two_grid/spectral_shape.hpp"

#include <algorithm>
#include <numeric>

#include <deal.II/lac/petsc_full_matrix.h>

namespace bart::calculator::two_grid {

namespace  {
using DealiiMatrix = dealii::FullMatrix<double>;
auto matrix_size_error_text(const DealiiMatrix& sigma_t, const DealiiMatrix& sigma_s) {
  return std::string{"Error in SpectralShape::CalculateSpectralShape, matrix size mismatch: Sigma_T matrix has "
                     "dimensions (" + std::to_string(sigma_t.m()) + ", " + std::to_string(sigma_t.n()) +
      "), and Sigma_S matrix has dimensions (" + std::to_string(sigma_s.m()) + ", "
                         + std::to_string(sigma_s.n())};
}
} // namespace

SpectralShape::SpectralShape(std::unique_ptr<EigenvalueSolver> eigenvalue_solver_ptr)
    : eigenvalue_solver_ptr_(std::move(eigenvalue_solver_ptr)) {
  AssertPointerNotNull(eigenvalue_solver_ptr_.get(), "eigenvalue solver", "SpectralShape constructor");
}

auto SpectralShape::CalculateSpectralShape(const DealiiMatrix& sigma_t,
                                           const DealiiMatrix& sigma_s) -> std::vector<double> {
  AssertThrow(sigma_t.m() == sigma_s.m() && sigma_t.n() == sigma_s.n() && sigma_t.m() == sigma_t.n() &&
      sigma_s.m() == sigma_s.n(), dealii::ExcMessage(matrix_size_error_text(sigma_t, sigma_s)))
  const int n_groups = sigma_t.m();
  DealiiMatrix downscattering(n_groups, n_groups), upscattering(n_groups, n_groups);

  for (int i = 0; i < n_groups; ++i) {
    for (int j = 0; j < i + 1; ++j)
      downscattering(i, j) = sigma_s(i, j);
    for (int j = i + 1; j < n_groups; ++j)
      upscattering(i, j) = sigma_s(i, j);
  }
  DealiiMatrix a(n_groups, n_groups);
  DealiiMatrix lhs(sigma_t);
  lhs.add(-1, downscattering);
  lhs.gauss_jordan();
  lhs.mmult(a, upscattering);

  dealii::PETScWrappers::FullMatrix petsc_matrix(n_groups, n_groups);
  for (int i = 0; i < n_groups; ++i) {
    for (int j = 0; j < n_groups; ++j) {
      petsc_matrix.set(i, j, a(i, j));
    }
  }
  petsc_matrix.compress(dealii::VectorOperation::insert);
  auto [eigenvalue, eigenvector] = this->eigenvalue_solver_ptr_->SpectralRadius(petsc_matrix);

  // Normalize in the L1 norm
  const double sum = std::accumulate(eigenvector.begin(), eigenvector.end(), 0.0,
                                     [](double running_sum, double val){ return running_sum + std::abs(val); });
  std::transform(eigenvector.begin(), eigenvector.end(), eigenvector.begin(),
                 [sum](const double val) { return std::abs(val) / sum; });

  return eigenvector;
}

} // namespace bart::calculator::two_grid
