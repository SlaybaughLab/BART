#ifndef BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
#define BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_

#include "formulation/types.h"
#include "formulation/equation/transport_scalar.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class Diffusion : public TransportScalar<dim> {
 public:
  Diffusion(const DiscretizationType discretization)
      : TransportScalar<dim>(discretization) {}

  virtual ~Diffusion() = default;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_