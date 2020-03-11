#include "framework/builder/framework_builder.h"

// Quadrature classes & factories
#include "quadrature/quadrature_generator_i.h"
#include "quadrature/factory/quadrature_factories.h"
#include "quadrature/utility/quadrature_utilities.h"

namespace bart {

namespace framework {

namespace builder {

template<int dim>
auto FrameworkBuilder<dim>::BuildQuadratureSet(
    const problem::ParametersI& problem_parameters)
-> std::shared_ptr<QuadratureSetType> {

  using QuadratureGeneratorType = quadrature::QuadratureGeneratorI<dim>;

  std::shared_ptr<QuadratureSetType> return_ptr = nullptr;
  std::shared_ptr<QuadratureGeneratorType > quadrature_generator_ptr = nullptr;

  const int order_value = problem_parameters.AngularQuadOrder();
  switch (problem_parameters.AngularQuad()) {
    default: {
      if (dim == 3) {
        quadrature_generator_ptr =
            quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
                quadrature::Order(order_value),
                quadrature::AngularQuadratureSetType::kLevelSymmetricGaussian);
      } else {
        AssertThrow(false,
                    dealii::ExcMessage("No supported quadratures for this dimension "
                                       "and transport model"))
      }
    }
  }

  return_ptr = quadrature::factory::MakeQuadratureSetPtr<dim>();

  auto quadrature_points = quadrature::utility::GenerateAllPositiveX<dim>(
      quadrature_generator_ptr->GenerateSet());

  quadrature::factory::FillQuadratureSet<dim>(return_ptr.get(),
                                              quadrature_points);

  return std::move(return_ptr);
}

template class FrameworkBuilder<1>;
template class FrameworkBuilder<2>;
template class FrameworkBuilder<3>;

} // namespace builder

} // namespace framework

} // namespace bart
