#include "framework/builder/cfem_framework_builder.h"

namespace bart {

namespace framework {

namespace builder {

template <int dim>
std::shared_ptr<formulation::CFEMStamperI> CFEM_FrameworkBuilder<dim>::BuildStamper(
    problem::ParametersI *problem_parameters,
    std::string material_mapping) {

}

template class CFEM_FrameworkBuilder<1>;
template class CFEM_FrameworkBuilder<2>;
template class CFEM_FrameworkBuilder<3>;

} // namespace builder

} // namespace framework

} // namespace bart