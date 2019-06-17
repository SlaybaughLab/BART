#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_

#include "quadrature/angular/angular_quadrature_set_i.h"

namespace bart {

namespace quadrature {

namespace angular {

template <int dim>
class AngularQuadratureSet : public AngularQuadratureSetI<dim> {
 public:
  AngularQuadratureSet() = default;
  ~AngularQuadratureSet() = default;

  std::vector<QuadraturePoint<dim>> quadrature_points() const override {
    return quadrature_points_;
  }

  int total_quadrature_points() const override {
    return quadrature_points_.size();
  }

 protected:
  std::vector<QuadraturePoint<dim>> quadrature_points_ = {};
};

} // namespace angular

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_