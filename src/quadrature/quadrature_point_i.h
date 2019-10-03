#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_

#include <memory>

#include "quadrature/ordinate_i.h"
#include "utility/named_type.h"

namespace bart {

namespace quadrature {

template <int dim>
class QuadraturePointI {
 public:
  virtual ~QuadraturePointI() = default;

  virtual std::shared_ptr<OrdinateI<dim>> ordinate() const = 0;
  virtual double weight() const = 0;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_
