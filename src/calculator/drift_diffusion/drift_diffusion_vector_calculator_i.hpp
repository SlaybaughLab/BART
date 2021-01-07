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
 * D(\phi, \vec{J}, \nabla\varphi, D) = \frac{1}{\phi}\left[\vec{J} - D\nabla\varphi\right]\;,
 * \f]
 *
 * where \f$\phi\f$ is the scalar flux, \f$J\f$ is the current, \f$\nabla\varphi\f$ is the gradient of the shape function,
 * and \f$D\f$ is the diffusion coefficent. The resulting return vector
 * will be of length \f$d \times 1\f$ where \f$d\f$ is the dimension of the problem.
 *
 * @tparam dim spatial dimension, required for return tensor and shape functions.
 *
 */
 template <int dim>
class DriftDiffusionVectorCalculatorI {
 public:
  using Tensor = typename dealii::Tensor<1, dim>;
  virtual auto DriftDiffusionVector(const double scalar_flux,
                                    const Tensor& current,
                                    const Tensor& shape_gradient,
                                    const double diffusion_coefficient) const -> Tensor = 0;
};

} // namespace bart::calculator::drift_diffusion

#endif //BART_SRC_CALCULATOR_DRIFT_DIFFUSION_DRIFT_DIFFUSION_VECTOR_CALCULATOR_I_HPP_
