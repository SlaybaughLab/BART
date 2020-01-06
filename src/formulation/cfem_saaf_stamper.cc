#include "formulation/cfem_saaf_stamper.h"

namespace bart {

namespace formulation {

template<int dim>
CFEM_SAAF_Stamper<dim>::CFEM_SAAF_Stamper(
    std::unique_ptr<SAAFFormulationType> saaf_ptr,
    std::shared_ptr<DomainDefinitionType> defintion_ptr)
    : formulation_ptr_(std::move(saaf_ptr)),
      definition_ptr_(defintion_ptr)
{}

template class CFEM_SAAF_Stamper<1>;
template class CFEM_SAAF_Stamper<2>;
template class CFEM_SAAF_Stamper<3>;

} // namespace formulation

} // namespace bart
