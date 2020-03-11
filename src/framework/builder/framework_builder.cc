#include "framework/builder/framework_builder.h"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>


// Convergence classes
#include "convergence/reporter/mpi_noisy.h"

// Domain classes
#include "domain/finite_element/finite_element_gaussian.h"

// Quadrature classes & factories
#include "quadrature/quadrature_generator_i.h"
#include "quadrature/factory/quadrature_factories.h"
#include "quadrature/utility/quadrature_utilities.h"

namespace bart {

namespace framework {

namespace builder {

template<int dim>
auto FrameworkBuilder<dim>::BuildConvergenceReporter()
-> std::unique_ptr<ReporterType> {
  using Reporter = bart::convergence::reporter::MpiNoisy;

  std::unique_ptr<ReporterType> return_ptr = nullptr;

  int this_process = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  auto pout_ptr = std::make_unique<dealii::ConditionalOStream>(std::cout, this_process == 0);
  return_ptr = std::make_unique<Reporter>(std::move(pout_ptr));

  return std::move(return_ptr);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildFiniteElement(ParametersType problem_parameters)
-> std::unique_ptr<FiniteElementType>{
  using FiniteElementGaussianType = domain::finite_element::FiniteElementGaussian<dim>;
  return std::make_unique<FiniteElementGaussianType>(
      problem::DiscretizationType::kContinuousFEM,
      problem_parameters.FEPolynomialDegree());
}

template<int dim>
auto FrameworkBuilder<dim>::BuildQuadratureSet(ParametersType problem_parameters)
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
