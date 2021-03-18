#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_I_HPP_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_I_HPP_

#include <memory>

#include "quadrature/ordinate_i.hpp"
#include "quadrature/quadrature_types.h"

namespace bart::quadrature {

/*! \brief Interface for quadrature points.
  * Classes derived from this class are designed to be quadrature points as part
 * of a quadrature set for collocation solving and numerical integration.
  * @tparam dim spatial dimension.
 */
template <int dim>
class QuadraturePointI {
 public:
  virtual ~QuadraturePointI() = default;

  /*! \brief Set the points ordinate and weight.
   * @return this object.
   */
  virtual auto SetTo(std::shared_ptr<OrdinateI<dim>>, quadrature::Weight) -> QuadraturePointI& = 0;

  /*! \brief Set the points ordinate.
   * @return this object.
   */
  virtual auto SetOrdinate(std::shared_ptr<OrdinateI<dim>>) -> QuadraturePointI& = 0;

  /*! \brief Set the points weight.
   * @return this object.
   */
  virtual auto SetWeight(quadrature::Weight) -> QuadraturePointI& = 0;

  /*! \brief Get the points ordinate.
   * @return a shared pointer to the ordinate.
   */
  virtual auto ordinate() const -> std::shared_ptr<OrdinateI<dim>> = 0;

  /*! \brief Get the points weight.
   * @return double value of the weight.
   */
  virtual auto weight() const -> double = 0;

  /*! \brief Get the cartesian position of the underlying ordinate.
   * @return array containing the cartesian position.
   */
  virtual auto cartesian_position() const -> std::array<double, dim> = 0;

  /*! \brief Get the cartesian position of the underlying ordinate as a tensor.
   * @return dealii tensor containing the caretisan position.
   */
   virtual auto cartesian_position_tensor() const -> dealii::Tensor<1, dim> = 0;
};

} // namespace bart::quadrature

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_I_HPP_
