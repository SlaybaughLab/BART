#ifndef BART_SRC_FORMULATION_CFEM_SAAF_STAMPER_H_
#define BART_SRC_FORMULATION_CFEM_SAAF_STAMPER_H_

#include "domain/definition_i.h"
#include "formulation/cfem_stamper_i.h"

#include "formulation/angular/cfem_self_adjoint_angular_flux_i.h"

#include <memory>

namespace bart {

namespace formulation {

template <int dim>
class CFEM_SAAF_Stamper {
 public:
  using DomainDefinitionType = domain::DefinitionI<dim>;
  using SAAFFormulationType = typename
      formulation::angular::CFEMSelfAdjointAngularFluxI<dim>;

   CFEM_SAAF_Stamper(std::unique_ptr<SAAFFormulationType> saaf_ptr,
                     std::shared_ptr<DomainDefinitionType> defintion_ptr);

  DomainDefinitionType* definition_ptr() const {return definition_ptr_.get();}
  SAAFFormulationType* formulation_ptr() const {return formulation_ptr_.get();}

 private:
  using InitializationTokenType =
  typename SAAFFormulationType::InitializationToken;

  // Dependencies
  std::unique_ptr<SAAFFormulationType> formulation_ptr_ = nullptr;
  std::shared_ptr<DomainDefinitionType> definition_ptr_ = nullptr;

  std::vector<formulation::CellPtr<dim>> cells_;
  InitializationTokenType saaf_initialization_token_;
};

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_CFEM_SAAF_STAMPER_H_
