#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_

#include <memory>

#include "quadrature/ordinate_i.h"
#include "quadrature/quadrature_types.h"

namespace bart {

namespace quadrature {

template <int dim>
class QuadraturePointI {
 public:
  virtual ~QuadraturePointI() = default;

  virtual QuadraturePointI& SetTo(
      const std::shared_ptr<OrdinateI<dim>>&, const quadrature::Weight) = 0;
  virtual QuadraturePointI& SetOrdinate(
      const std::shared_ptr<OrdinateI<dim>>&) = 0;
  virtual QuadraturePointI& SetWeight(const quadrature::Weight) = 0;

  virtual std::shared_ptr<OrdinateI<dim>> ordinate() const = 0;
  virtual double weight() const = 0;
  virtual std::array<double, dim> cartesian_position() const = 0;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_I_H_
