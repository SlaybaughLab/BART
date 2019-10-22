#ifndef BART_SRC_QUADRATURE_ORDINATE_I_H_
#define BART_SRC_QUADRATURE_ORDINATE_I_H_

#include <deal.II/base/tensor.h>

#include "quadrature/quadrature_types.h"
#include "utility/named_type.h"

namespace bart {

namespace quadrature {

/*! \brief Interface for ordinates for quadrature sets.
 *
 * @tparam dim spatial dimension.
 */
template <int dim>
class OrdinateI {
 public:
  virtual ~OrdinateI() = default;
  virtual OrdinateI& set_cartesian_position(const CartesianPosition<dim>) = 0;
  virtual std::array<double, dim> cartesian_position() const = 0;
  virtual dealii::Tensor<1, dim> cartesian_position_tensor() const = 0;

  virtual bool operator==(const OrdinateI&) const = 0;
  virtual bool operator!=(const OrdinateI&) const = 0;
  virtual bool operator==(const std::array<double, dim>) const = 0;
  virtual bool operator!=(const std::array<double, dim>) const = 0;
};

template<int dim>
std::array<double, dim> Reflect(const OrdinateI<dim>& ordinate) {
  std::array<double, dim> return_array;
  for (unsigned int i = 0; i < dim; ++i)
    return_array.at(i) = -ordinate.cartesian_position().at(i);
  return return_array;
}

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_ORDINATE_I_H_
