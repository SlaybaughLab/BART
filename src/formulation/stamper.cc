#include "formulation/stamper.h"

namespace bart {

namespace formulation {

template<int dim>
Stamper<dim>::Stamper(std::shared_ptr<domain::DefinitionI<dim>>)
    : domain_ptr_(nullptr) {}

template class Stamper<1>;
template class Stamper<2>;
template class Stamper<3>;

} // namespace formulation

} // namespace bart
