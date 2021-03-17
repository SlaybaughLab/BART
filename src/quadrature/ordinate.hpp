#ifndef BART_SRC_QUADRATURE_ORDINATE_HPP_
#define BART_SRC_QUADRATURE_ORDINATE_HPP_

#include "quadrature/ordinate_i.hpp"
#include "quadrature/quadrature_types.h"

namespace bart::quadrature {

/*! \brief Default implementation of Ordinate class.
 *
 * @tparam dim  spatial dimension.
 */
template <int dim>
class Ordinate : public OrdinateI<dim> {
 public:
  /*! \brief Default constructor, cartesian position is the origin. */
  Ordinate() { cartesian_position_.fill(0); }
  /*! \brief Constructor based on provided cartesian position. */
  explicit Ordinate(const CartesianPosition<dim> position) : cartesian_position_(std::move(position.get())) {}
  virtual ~Ordinate() = default;

  auto set_cartesian_position(const CartesianPosition<dim> to_set) -> Ordinate& override {
    cartesian_position_ = to_set.get();
    return *this;
  }

  auto cartesian_position() const -> std::array<double, dim> override { return cartesian_position_; }
  auto cartesian_position_tensor() const -> dealii::Tensor<1, dim> override {
    dealii::Tensor<1, dim> tensor_position;
    for (int i = 0; i < dim; ++i)
      tensor_position[i] = cartesian_position_.at(i);
    return tensor_position;
  }

  auto operator==(const OrdinateI<dim>& rhs) const -> bool override {
    return cartesian_position_ == rhs.cartesian_position();
  }
  auto operator!=(const OrdinateI<dim>& rhs) const -> bool override { return !operator==(rhs); }
  auto operator==(const std::array<double, dim> rhs) const -> bool override { return cartesian_position_ == rhs; };
  auto operator!=(const std::array<double, dim> rhs) const -> bool override { return cartesian_position_ != rhs; };

 private:
  std::array<double, dim> cartesian_position_;
};

} // namespace bart::quadrature

#endif //BART_SRC_QUADRATURE_ORDINATE_HPP_
