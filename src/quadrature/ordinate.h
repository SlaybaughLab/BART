#ifndef BART_SRC_QUADRATURE_ORDINATE_H_
#define BART_SRC_QUADRATURE_ORDINATE_H_

#include "quadrature/ordinate_i.h"

namespace bart {

namespace quadrature {

template <int dim>
using CartesianPosition = utility::NamedType<std::array<double, dim>,
                                             struct CartesianPositionParameter>;

/*! \brief Default implementation of Ordinate class.
 *
 * @tparam dim  spatial dimension.
 */
template <int dim>
class Ordinate : public OrdinateI<dim> {
 public:
  explicit Ordinate(CartesianPosition<dim>);
  virtual ~Ordinate() = default;

  std::array<double, dim> cartesian_position() const override {
    return cartesian_position_;
  }

  dealii::Tensor<1, dim> cartesian_position_tensor() const override {
    dealii::Tensor<1, dim> return_tensor;
    for (int i = 0; i < dim; ++i)
      return_tensor[i] = cartesian_position_.at(i);
    return return_tensor;
  }



 private:
  std::array<double, dim> cartesian_position_;
};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_ORDINATE_H_
