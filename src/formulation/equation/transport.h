#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_

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

  EquationType equation_type() const override { return equation_type_; };
  DiscretizationType discretization_type() const override {
    return discretization_type_; };

 protected:
  EquationType equation_type_;
  DiscretizationType discretization_type_;

};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_