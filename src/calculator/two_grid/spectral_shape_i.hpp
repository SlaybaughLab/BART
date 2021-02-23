#ifndef BART_SRC_CALCULATOR_TWO_GRID_SPECTRAL_SHAPE_I_HPP_
#define BART_SRC_CALCULATOR_TWO_GRID_SPECTRAL_SHAPE_I_HPP_

#include <deal.II/lac/full_matrix.h>

#include "utility/has_description.h"

//! Calculators for the two-grid acceleration method <a href="https://doi.org/10.13182/NSE115-253">Adams and Morel (2017)</a>
namespace bart::calculator::two_grid {

/*! \brief Interface for classes that calculate the spectral shape function for the two-grid acceleration method.
 *
 * The spectral shape function is defined for the two-grid method of <a href="https://doi.org/10.13182/NSE115-253">Adams and Morel (2017)</a>.
 * This function must be calculated in each homogenous material region, and is calculated by solving
 * for the eigenvector \f$\xi \in \mathbb{R}^{(G \times 1)}\f$ in the following equation:
 * \f[
 * \left(\bf{T} - \bf{S}_{D0}\right)^{-1}\bf{S}_{U0}\vec{\xi} = \rho\vec{\xi}\;,
 * \f]
 * where \f$\bf{T}\f$ is the total cross-section matrix, \f$\bf{S}_{D0}\f$ is the isotropic downscatter matrix,
 * \f$\bf{S}_{U0}\f$ is the isotropic upscatter matrix, and \f$\rho\f$ is the spectral radius. The spectral shape
 * function is normalized such that
 * \f[
 * \sum_{g = 0}^{G - 1}\xi_g = 1\;.
 * \f]
 *
 */
class SpectralShapeI : public utility::HasDescription {
 public:
  using DealiiMatrix = dealii::FullMatrix<double>;
  virtual ~SpectralShapeI() = default;
  /*! \brief Calculate the spectral shape.
   *
   * @param sigma_t total cross-section matrix.
   * @param sigma_s full scattering matrix.
   * @return spectral shape function
   */
  virtual auto CalculateSpectralShape(const DealiiMatrix& sigma_t,
                                      const DealiiMatrix& sigma_s) -> std::vector<double> = 0;
};

} // namespace bart::calculator::two_grid

#endif //BART_SRC_CALCULATOR_TWO_GRID_SPECTRAL_SHAPE_I_HPP_
