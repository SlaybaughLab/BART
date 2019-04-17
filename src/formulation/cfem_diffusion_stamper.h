#ifndef BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_
#define BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_

#include <memory>
#include <vector>

#include "domain/definition_i.h"
#include "formulation/scalar/cfem_diffusion_i.h"
#include "formulation/stamper_i.h"

namespace bart {

namespace formulation {

template <int dim>
class CFEM_DiffusionStamper : public StamperI<dim> {
 public:
  CFEM_DiffusionStamper(
      std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
      std::unique_ptr<domain::DefinitionI<dim>> definition_ptr);

 private:
  using InitializationToken =
      typename formulation::scalar::CFEM_DiffusionI<dim>::InitializationToken;
  using Cell = typename domain::DefinitionI<dim>::Cell;

  std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr_;
  std::unique_ptr<domain::DefinitionI<dim>> definition_ptr_;
  InitializationToken diffusion_init_token_;
  std::vector<Cell> cells_;
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_CFEM_DIFFUSION_STAMPER_H_