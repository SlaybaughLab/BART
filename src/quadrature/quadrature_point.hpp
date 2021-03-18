#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_HPP_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_HPP_

#include "quadrature/quadrature_point_i.hpp"
#include "quadrature/quadrature_types.h"

namespace bart::quadrature {

/*! \brief Default implementation of the QuadraturePoint interface.
 * @tparam dim spatial dimension.
 */
template <int dim>
class QuadraturePoint : public QuadraturePointI<dim> {
 public:
  /// \brief Default constructor.
  QuadraturePoint() = default;

  /*! \brief Constructor with provided ordinate and weight.
   * @param ordinate the quadrature point's ordinate.
   * @param weight the quadrature point's weight.
   */
  QuadraturePoint(std::shared_ptr<OrdinateI<dim>> ordinate, Weight weight);

  auto SetTo(std::shared_ptr<OrdinateI<dim>>, Weight) -> QuadraturePoint<dim>& override;
  auto SetOrdinate(std::shared_ptr<OrdinateI<dim>>) -> QuadraturePoint<dim>& override;
  auto SetWeight(Weight) -> QuadraturePoint<dim>& override;

  [[nodiscard]] auto ordinate() const -> std::shared_ptr<OrdinateI<dim>> override {return ordinate_;}
  [[nodiscard]] auto weight() const -> double override { return weight_; }
  [[nodiscard]] auto cartesian_position() const -> std::array<double, dim> override {
    return ordinate_->cartesian_position(); }
  [[nodiscard]] auto cartesian_position_tensor() const -> dealii::Tensor<1, dim> override {
    return ordinate_->cartesian_position_tensor(); }
 private:
  std::shared_ptr<OrdinateI<dim>> ordinate_{ nullptr };
  double weight_{ 0 };
};

} // namespace bart::quadrature

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_HPP_
