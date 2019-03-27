#ifndef BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_H_
#define BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_H_

#include <memory>

#include "data/cross_sections.h"
#include "domain/finite_element_i.h"
#include "formulation/scalar/cfem_i.h"

namespace bart {

namespace formulation {

namespace scalar {

template <int dim>
class CFEM_Diffusion : public CFEM_I {
 public:
  CFEM_Diffusion(std::shared_ptr<domain::FiniteElementI<dim>> finite_element,
                 std::shared_ptr<data::CrossSections> cross_sections);

 protected:
};


} // namespace scalar

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_H_