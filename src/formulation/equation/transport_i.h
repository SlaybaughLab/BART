#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_I_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_I_H_

#include <memory>

#include "formulation/types.h"
#include "domain/finite_element_i.h"


namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class TransportI {
 public:
  virtual ~TransportI() = default;

  virtual EquationType equation_type() const = 0;
  virtual DiscretizationType discretization_type() const = 0;
  virtual TransportI& ProvideFiniteElement(
      std::shared_ptr<domain::FiniteElementI<dim>>);
};

} // namespace formulation

} // namespace equation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_