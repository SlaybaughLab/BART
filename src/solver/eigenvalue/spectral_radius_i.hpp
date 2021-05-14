#ifndef BART_SRC_SOLVER_EIGENVALUE_SPECTRAL_RADIUS_I_HPP_
#define BART_SRC_SOLVER_EIGENVALUE_SPECTRAL_RADIUS_I_HPP_

#include <deal.II/lac/petsc_matrix_base.h>

namespace bart::solver::eigenvalue {

class SpectralRadiusI {
 public:
  using MatrixBase = dealii::PETScWrappers::MatrixBase;
  virtual ~SpectralRadiusI() = default;
  /*! \brief Calculate the spectral radius and corresponding eigenvector.
   *
   * @return pair containing the spectral radius and eigenvector.
   */
  virtual auto SpectralRadius(const MatrixBase&) -> std::pair<double, std::vector<double>> = 0;
};

} // namespace bart::solver::eigenvalue

#endif //BART_SRC_SOLVER_EIGENVALUE_SPECTRAL_RADIUS_I_HPP_
