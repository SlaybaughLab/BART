#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_FINITE_ELEMENT_FACTORY_H_
#define BART_SRC_DOMAIN_FINITE_ELEMENT_FINITE_ELEMENT_FACTORY_H_

#include <memory>

#include "problem/parameter_types.h"

namespace bart {

namespace domain {

namespace finite_element {
template <int dim> class FiniteElementI;

enum class FiniteElementImpl {
  kGaussian = 0,
};

template <int dim>
class FiniteElementFactory {
 public:
  static std::unique_ptr<FiniteElementI<dim>> MakeFiniteElement(
      const problem::DiscretizationType discretization_type,
      const int polynomial_degree,
      const FiniteElementImpl finite_element_type);
};

} // namespace finite_element

} // namespace domain

} //namespace bart

#endif //BART_SRC_DOMAIN_FINITE_ELEMENT_FINITE_ELEMENT_FACTORY_H_
