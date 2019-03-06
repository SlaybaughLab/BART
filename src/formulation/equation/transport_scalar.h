#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_

#include "formulation/types.h"
#include "formulation/equation/transport.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class TransportScalar : public Transport<dim> {
 public:
  virtual ~TransportScalar() = default;

  TransportScalar(const DiscretizationType discretization)
      : Transport<dim>(EquationType::kScalar, discretization) {};

  ScalarEquations scalar_equation() const { return scalar_equation_; };
 protected:
  ScalarEquations scalar_equation_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_