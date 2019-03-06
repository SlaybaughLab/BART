#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_

#include <memory>

#include "domain/finite_element_i.h"
#include "formulation/types.h"
#include "formulation/equation/transport_i.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class Transport : public TransportI<dim> {
 public:
  Transport(EquationType equation_type, DiscretizationType discretization_type)
      : equation_type_(equation_type),
        discretization_type_(discretization_type) {}
  virtual ~Transport() = default;

  Transport& ProvideFiniteElement(
      std::shared_ptr<domain::FiniteElementI<dim>> finite_element) override {
    finite_element_ = finite_element;
  }

  EquationType equation_type() const override { return equation_type_; };
  DiscretizationType discretization_type() const override {
    return discretization_type_; };

 protected:
  EquationType equation_type_;
  DiscretizationType discretization_type_;
  std::shared_ptr<domain::FiniteElementI<dim>> finite_element_;

};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_