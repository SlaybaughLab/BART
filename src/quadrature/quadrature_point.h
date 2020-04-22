#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_H_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_H_

#include "quadrature/quadrature_point_i.h"
#include "quadrature/quadrature_types.h"

namespace bart {

namespace quadrature {

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

  QuadraturePoint<dim> &SetTo(const std::shared_ptr<OrdinateI<dim>> &,
                              const Weight) override;
  QuadraturePoint<dim> &SetOrdinate(
      const std::shared_ptr<OrdinateI<dim>> &) override;
  QuadraturePoint<dim> &SetWeight(const Weight) override;

  std::shared_ptr<OrdinateI<dim>> ordinate() const override {return ordinate_;}
  double weight() const override { return weight_; }
  std::array<double, dim> cartesian_position() const override {
    return ordinate_->cartesian_position(); }
  dealii::Tensor<1, dim> cartesian_position_tensor() const override {
    return ordinate_->cartesian_position_tensor(); }
 private:
  std::shared_ptr<OrdinateI<dim>> ordinate_ = nullptr;
  double weight_ = 0;
};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_H_
