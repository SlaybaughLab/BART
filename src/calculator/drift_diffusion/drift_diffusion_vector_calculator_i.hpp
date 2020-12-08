#ifndef BART_SRC_CALCULATOR_DRIFT_DIFFUSION_DRIFT_DIFFUSION_VECTOR_CALCULATOR_I_HPP_
#define BART_SRC_CALCULATOR_DRIFT_DIFFUSION_DRIFT_DIFFUSION_VECTOR_CALCULATOR_I_HPP_

#include <deal.II/base/tensor.h>

namespace bart::calculator::drift_diffusion {

/*! \brief Interface for classes that calculate the drift-diffusion vector.
 *
 * This class has the single job of combining the terms of the drift-diffusion vector and returning the correct value.
 * The vector is calculated for a single degree of freedom and quadrature point within the cell. Therefore, the values
 * of the scalar flux and integrated angular flux are both scalar values.
 *
 * \f[
 * D(\phi, I(\Psi), \nabla\varphi, \sigma_t, D) = \frac{1}{\phi}\left[\frac{1}{\sigma_t}I(\Psi)\nabla\varphi
 *                                                - D\nabla\varphi\phi\right]\;,
 * \f]
 *
 * where \f$\phi\f$ is the scalar flux, \f$I(\Psi)\f$ is the angular flux integrated over angle (using quadrature)
 * \f$\sum_m w_m\hat{\Omega}_m\hat{\Omega_m}\Psi_m\f$, \f$\nabla\varphi\f$ is the gradient of the shape function,
 * \f$\sigma_t\f$ is the total cross-section, and \f$D\f$ is the diffusion coefficent. The resulting return vector
 * will be of length \f$d \times 1\f$ where \f$d\f$ is the dimension of the problem.
 *
 * @tparam dim spatial dimension, required for return tensor and shape functions.
 *
 */
 template <int dim>
class DriftDiffusionVectorCalculatorI {
 public:
  using Tensor = typename dealii::Tensor<1, dim>;
  virtual auto DriftDiffusion(const double scalar_flux,
                              const double integrated_angular_flux,
                              const Tensor& shape_gradient,
                              const double sigma_t,
                              const double diffusion_coefficient) const -> Tensor = 0;
};

} // namespace bart::calculator::drift_diffusion

#endif //BART_SRC_CALCULATOR_DRIFT_DIFFUSION_DRIFT_DIFFUSION_VECTOR_CALCULATOR_I_HPP_
