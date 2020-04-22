#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_

#include <memory>

#include "quadrature/ordinate_i.h"
#include "quadrature/quadrature_types.h"

namespace bart {

namespace quadrature {

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
  virtual QuadraturePointI& SetTo(
      const std::shared_ptr<OrdinateI<dim>>&, const quadrature::Weight) = 0;

  /*! \brief Set the points ordinate.
   * @return this object.
   */
  virtual QuadraturePointI& SetOrdinate(
      const std::shared_ptr<OrdinateI<dim>>&) = 0;

  /*! \brief Set the points weight.
   * @return this object.
   */
  virtual QuadraturePointI& SetWeight(const quadrature::Weight) = 0;

  /*! \brief Get the points ordinate.
   * @return a shared pointer to the ordinate.
   */
  virtual std::shared_ptr<OrdinateI<dim>> ordinate() const = 0;

  /*! \brief Get the points weight.
   * @return double value of the weight.
   */
  virtual double weight() const = 0;

  /*! \brief Get the cartesian position of the underlying ordinate.
   * @return array containing the cartesian position.
   */
  virtual std::array<double, dim> cartesian_position() const = 0;

  /*! \brief Get the cartesian position of the underlying ordinate as a tensor.
   * @return dealii tensor containing the caretisan position.
   */
   virtual dealii::Tensor<1, dim> cartesian_position_tensor() const = 0;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_
