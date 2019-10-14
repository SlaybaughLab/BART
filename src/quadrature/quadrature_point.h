#ifndef BART_SRC_QUADRATURE_QUADRATURE_POINT_H_
#define BART_SRC_QUADRATURE_QUADRATURE_POINT_H_

#include "quadrature/quadrature_point_i.h"
#include "quadrature/quadrature_types.h"

namespace bart {

namespace quadrature {

template <int dim>
class QuadraturePoint : public QuadraturePointI<dim> {
 public:
  QuadraturePoint() = default;
  QuadraturePoint(std::shared_ptr<OrdinateI<dim>> ordinate, Weight);

  QuadraturePoint<dim> &SetTo(const std::shared_ptr<OrdinateI<dim>> &,
                              const Weight) override;
  QuadraturePoint<dim> &SetOrdinate(
      const std::shared_ptr<OrdinateI<dim>> &) override;
  QuadraturePoint<dim> &SetWeight(const Weight) override;

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
