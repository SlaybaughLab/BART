#include "framework/builder/framework_builder.h"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <sstream>

#include "utility/reporter/mpi.h"
#include "utility/reporter/colors.h"

// Convergence classes
#include "convergence/final_checker_or_n.h"
#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/parameters/single_parameter_checker.h"
#include "convergence/reporter/mpi.h"

// Domain classes
#include "domain/definition.h"
#include "domain/finite_element/finite_element_gaussian.h"
#include "domain/mesh/mesh_cartesian.h"

// Formulation classes
#include "formulation/angular/self_adjoint_angular_flux.h"
#include "formulation/scalar/diffusion.h"
#include "formulation/stamper.h"
#include "formulation/updater/saaf_updater.h"
#include "formulation/updater/diffusion_updater.h"

// Material classes
#include "material/material_protobuf.h"

// Solver classes
#include "solver/group/single_group_solver.h"
#include "solver/gmres.h"

// Iteration classes
#include "iteration/initializer/initialize_fixed_terms_once.h"

// Quadrature classes & factories
#include "quadrature/quadrature_generator_i.h"
#include "quadrature/factory/quadrature_factories.h"
#include "quadrature/utility/quadrature_utilities.h"

namespace bart {

namespace framework {

namespace builder {

template<int dim>
void FrameworkBuilder<dim>::BuildFramework(std::string name,
                                           ParametersType& prm) {

  reporter_ptr_->Report("Building Framework: ");
  reporter_ptr_->Report(name + "\n", utility::reporter::Color::Green);
  std::shared_ptr<FiniteElementType> finite_element_ptr =
      std::move(BuildFiniteElement(prm));
  std::shared_ptr<CrossSectionType> cross_sections_ptr =
      std::move(BuildCrossSections(prm));

  auto domain_ptr = BuildDomain(prm, finite_element_ptr, "");

}

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
auto FrameworkBuilder<dim>::BuildCrossSections(
    const problem::ParametersI& problem_parameters)
    -> std::unique_ptr<CrossSectionType> {
  reporter_ptr_->Report("\tBuilding Cross-sections: ");
  std::unique_ptr<CrossSectionType> return_ptr = nullptr;
  // Default implementation using protocol buffers
  try {
    MaterialProtobuf materials(problem_parameters.MaterialFilenames(),
                               problem_parameters.IsEigenvalueProblem(),
                               problem_parameters.DoNDA(),
                               problem_parameters.NEnergyGroups(),
                               problem_parameters.NumberOfMaterials());
    return_ptr = std::move(std::make_unique<CrossSectionType>(materials));
    reporter_ptr_->Report("Built (default) Cross-sections using protobuf\n",
                          utility::reporter::Color::Green);
  } catch (...) {
    reporter_ptr_->Report("Error building (default) Cross-sections using protobuf\n",
                          utility::reporter::Color::Red);
    throw;
  }
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildDiffusionFormulation(
    const std::shared_ptr<FiniteElementType>& finite_element_ptr,
    const std::shared_ptr<data::CrossSections>& cross_sections_ptr,
    const formulation::DiffusionFormulationImpl implementation)
-> std::unique_ptr<DiffusionFormulationType> {
  std::unique_ptr<DiffusionFormulationType> return_ptr = nullptr;

  if (implementation == formulation::DiffusionFormulationImpl::kDefault) {
    using ReturnType = formulation::scalar::Diffusion<dim>;
    return_ptr = std::move(std::make_unique<ReturnType>(
        finite_element_ptr, cross_sections_ptr));
  }

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildDomain(
    ParametersType problem_parameters,
    const std::shared_ptr<FiniteElementType>& finite_element_ptr,
    std::string material_mapping)
-> std::unique_ptr<DomainType>{
  reporter_ptr_->Report("\tBuilding Mesh: ");
  auto mesh_ptr = std::make_unique<domain::mesh::MeshCartesian<dim>>(
      problem_parameters.SpatialMax(),
      problem_parameters.NCells(),
      material_mapping);
  reporter_ptr_->Report("Built: " + mesh_ptr->description() + "\n",
                        utility::reporter::Color::Green);

  return std::make_unique<domain::Definition<dim>>(std::move(mesh_ptr),
                                                   finite_element_ptr);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildFiniteElement(ParametersType problem_parameters)
-> std::unique_ptr<FiniteElementType>{
  std::unique_ptr<FiniteElementType> return_ptr = nullptr;

  using FiniteElementGaussianType = domain::finite_element::FiniteElementGaussian<dim>;

  reporter_ptr_->Report("\tBuilding Finite Element: ");

  try {
    return_ptr = std::move(std::make_unique<FiniteElementGaussianType>(
        problem::DiscretizationType::kContinuousFEM,
        problem_parameters.FEPolynomialDegree()));

    reporter_ptr_->Report("Built: " + return_ptr->description() + "\n",
                          utility::reporter::Color::Green);
  } catch (...) {
    reporter_ptr_->Report("Error");
    throw;
  }
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildFixedUpdater(
    std::unique_ptr<DiffusionFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr)
-> std::unique_ptr<FixedUpdaterType> {
  std::unique_ptr<FixedUpdaterType> return_ptr = nullptr;

  using ReturnType = formulation::updater::DiffusionUpdater<dim>;
  return_ptr = std::move(std::make_unique<ReturnType>(
      std::move(formulation_ptr),
      std::move(stamper_ptr)));

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildFixedUpdater(
    std::unique_ptr<SAAFFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr,
    const std::shared_ptr<QuadratureSetType>& quadrature_set_ptr)
-> std::unique_ptr<FixedUpdaterType> {
  std::unique_ptr<FixedUpdaterType> return_ptr = nullptr;

  using ReturnType = formulation::updater::SAAFUpdater<dim>;
  return_ptr = std::move(std::make_unique<ReturnType>(
      std::move(formulation_ptr),
      std::move(stamper_ptr),
      quadrature_set_ptr));

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildInitializer(
    const std::shared_ptr<formulation::updater::FixedUpdaterI>& updater_ptr,
    const int total_groups,
    const int total_angles) -> std::unique_ptr<InitializerType> {
  std::unique_ptr<InitializerType> return_ptr = nullptr;

  using InitializeOnceType = iteration::initializer::InitializeFixedTermsOnce;

  return_ptr = std::move(std::make_unique<InitializeOnceType>(
      updater_ptr, total_groups, total_angles));

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildMomentConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<MomentConvergenceCheckerType>{
  //TODO(Josh): Add option for using other than L1Norm

  using CheckerType = convergence::moments::SingleMomentCheckerL1Norm;
  using FinalCheckerType = convergence::FinalCheckerOrN<
      system::moments::MomentVector,
      convergence::moments::SingleMomentCheckerI>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(single_checker_ptr));

  return_ptr->SetMaxIterations(max_iterations);

  return std::move(return_ptr);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildParameterConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<ParameterConvergenceCheckerType>{

  using CheckerType = convergence::parameters::SingleParameterChecker;
  using FinalCheckerType = convergence::FinalCheckerOrN<double, CheckerType>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(single_checker_ptr));

  return_ptr->SetMaxIterations(max_iterations);

  return std::move(return_ptr);
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

template <int dim>
auto FrameworkBuilder<dim>::BuildSAAFFormulation(
    const std::shared_ptr<FiniteElementType>& finite_element_ptr,
    const std::shared_ptr<data::CrossSections>& cross_sections_ptr,
    const std::shared_ptr<QuadratureSetType>& quadrature_set_ptr,
    const formulation::SAAFFormulationImpl implementation)
-> std::unique_ptr<SAAFFormulationType> {
  std::unique_ptr<SAAFFormulationType> return_ptr;

  if (implementation == formulation::SAAFFormulationImpl::kDefault) {
    using ReturnType = formulation::angular::SelfAdjointAngularFlux<dim>;

    return_ptr = std::move(std::make_unique<ReturnType>(finite_element_ptr,
                                                        cross_sections_ptr,
                                                        quadrature_set_ptr));
  }

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildSingleGroupSolver(
    const int max_iterations,
    const double convergence_tolerance)
-> std::unique_ptr<SingleGroupSolverType> {
  std::unique_ptr<SingleGroupSolverType> return_ptr = nullptr;

  auto linear_solver_ptr = std::make_unique<solver::GMRES>(max_iterations,
                                                           convergence_tolerance);

  return_ptr = std::move(std::make_unique<solver::group::SingleGroupSolver>(
          std::move(linear_solver_ptr)));

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildStamper(
    const std::shared_ptr<DomainType>& domain_ptr)
-> std::unique_ptr<StamperType> {
  std::unique_ptr<StamperType> return_ptr = nullptr;

  return_ptr = std::move(
      std::make_unique<formulation::Stamper<dim>>(domain_ptr));

  return return_ptr;
}

template class FrameworkBuilder<1>;
template class FrameworkBuilder<2>;
template class FrameworkBuilder<3>;

} // namespace builder

} // namespace framework

} // namespace bart
