#include "formulation/cfem_diffusion_stamper.h"

namespace bart {

namespace formulation {

template<int dim>
CFEM_DiffusionStamper<dim>::CFEM_DiffusionStamper(
    std::unique_ptr<formulation::scalar::CFEM_DiffusionI<dim>> diffusion_ptr,
    std::unique_ptr<domain::DefinitionI<dim>> domain_ptr)
    : diffusion_ptr_(std::move(diffusion_ptr)),
      domain_ptr_(std::move(domain_ptr)) {


      }

template class CFEM_DiffusionStamper<1>;
template class CFEM_DiffusionStamper<2>;
template class CFEM_DiffusionStamper<3>;

} // namespace formulation

} // namespace bart