#include "formulation/cfem_diffusion_stamper.h"

namespace bart {

namespace formulation {

template<int dim>
CFEM_DiffusionStamper<dim>::CFEM_DiffusionStamper(
    std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
    std::unique_ptr<domain::DefinitionI<dim>> definition_ptr)
    : diffusion_ptr_(std::move(diffusion_ptr)),
      definition_ptr_(std::move(definition_ptr)) {

      cells_ = definition_ptr_->Cells();
      diffusion_init_token_ = diffusion_ptr_->Precalculate(cells_[0]);
}

template class CFEM_DiffusionStamper<1>;
template class CFEM_DiffusionStamper<2>;
template class CFEM_DiffusionStamper<3>;

} // namespace formulation

} // namespace bart