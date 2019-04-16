#ifndef BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_
#define BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_

#include <memory>

#include "domain/definition_i.h"
#include "formulation/stamper_i.h"

namespace bart {

namespace formulation {

template <int dim>
class CFEM_DiffusionStamper : public StamperI<dim> {
 public:
  explicit CFEM_DiffusionStamper(std::unique_ptr<domain::DefinitionI<dim>> domain_ptr)
      : domain_ptr_(std::move(domain_ptr)) {}

 private:
  std::unique_ptr<domain::DefinitionI<dim>> domain_ptr_;
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_