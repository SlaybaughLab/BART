#include "framework/builder/cfem_framework_builder.h"

#include "problem/parameter_types.h"
#include "domain/finite_element_gaussian.h"
#include "domain/mesh_cartesian.h"

namespace bart {

namespace framework {

namespace builder {

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildFiniteElement(
    problem::ParametersI *problem_parameters)-> std::shared_ptr<FiniteElement> {
  return std::make_shared<domain::FiniteElementGaussian<dim>>(
      problem::DiscretizationType::kContinuousFEM,
      problem_parameters->FEPolynomialDegree());
}

template <int dim>
std::shared_ptr<formulation::CFEMStamperI> CFEM_FrameworkBuilder<dim>::BuildStamper(
    problem::ParametersI *problem_parameters,
    std::string material_mapping) {

  auto mesh_ptr = std::make_unique<domain::MeshCartesian<dim>>(
      problem_parameters->SpatialMax(),
      problem_parameters->NCells(),
      material_mapping
      );

  auto finite_element_ptr = std::make_unique<domain::FiniteElementGaussian<dim>>(
      problem::DiscretizationType::kContinuousFEM,
      problem_parameters->FEPolynomialDegree()
      );

}


template class CFEM_FrameworkBuilder<1>;
template class CFEM_FrameworkBuilder<2>;
template class CFEM_FrameworkBuilder<3>;

} // namespace builder

} // namespace framework

} // namespace bart