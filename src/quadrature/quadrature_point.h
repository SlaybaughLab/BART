#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_H_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_H_

#include "quadrature/quadrature_point_i.h"

namespace bart {

namespace quadrature {

template <int dim>
class QuadraturePoint : public QuadraturePointI<dim> {
 public:

  QuadraturePoint(std::shared_ptr<OrdinateI<dim>> ordinate, double weight);

  std::shared_ptr<OrdinateI<dim>> ordinate() const override {
    return ordinate_;
  }
  double weight() const override {
    return weight_;
  }
 private:
  std::shared_ptr<OrdinateI<dim>> ordinate_ = nullptr;
  double weight_ = 0;
};

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_POINT_H_
