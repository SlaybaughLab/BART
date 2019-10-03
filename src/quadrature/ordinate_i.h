#ifndef BART_SRC_QUADRATURE_ORDINATE_I_H_
#define BART_SRC_QUADRATURE_ORDINATE_I_H_

#include <deal.II/base/tensor.h>

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
  virtual std::array<double, dim> cartesian_position() const = 0;
  virtual dealii::Tensor<1, dim> cartesian_position_tensor() const = 0;

  virtual bool operator==(const OrdinateI&) const = 0;
  virtual bool operator!=(const OrdinateI&) const = 0;
  virtual bool operator==(const std::array<double, dim>) const = 0;
  virtual bool operator!=(const std::array<double, dim>) const = 0;
  //virtual OrdinateI<dim> operator-() const = 0;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_ORDINATE_I_H_
