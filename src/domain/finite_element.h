#ifndef BART_DOMAIN_FINITE_ELEMENT_H_
#define BART_DOMAIN_FINITE_ELEMENT_H_

#include "domain/finite_element_i.h"

namespace bart {

namespace domain {

template<int dim>
class FiniteElement : public FiniteElementI<dim> {
 public:
  using typename FiniteElementI<dim>::CellPtr;
  virtual ~FiniteElement() = default;

  bool SetCell(const CellPtr &to_set) override {

  }

};


} // namespace domain

} // namespace bart

#endif // BART_DOMAIN_FINITE_ELEMENT_H_