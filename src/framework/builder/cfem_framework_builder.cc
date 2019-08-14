#include "framework/builder/cfem_framework_builder.h"

#include "domain/definition.h"
#include "domain/finite_element/finite_element_gaussian.h"
#include "domain/mesh/mesh_cartesian.h"
#include "formulation/cfem_diffusion_stamper.h"
#include "formulation/scalar/cfem_diffusion.h"
#include "iteration/updater/source_updater_gauss_seidel.h"
#include "iteration/updater/fixed_updater.h"
#include "problem/parameter_types.h"
#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/parameters/single_parameter_checker.h"
#include "convergence/final_checker_or_n.h"


namespace bart {

namespace framework {

namespace builder {

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildFiniteElement(
    problem::ParametersI *problem_parameters)-> std::unique_ptr<FiniteElement> {
  return std::make_unique<domain::finite_element::FiniteElementGaussian<dim>>(
      problem::DiscretizationType::kContinuousFEM,
      problem_parameters->FEPolynomialDegree());
}

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildDomain(
    problem::ParametersI *problem_parameters,
    const std::shared_ptr<FiniteElement> &finite_element_ptr,
    std::string material_mapping)-> std::unique_ptr<Domain> {

  // Build mesh
  auto mesh_ptr = std::make_unique<domain::mesh::MeshCartesian<dim>>(
      problem_parameters->SpatialMax(),
          problem_parameters->NCells(),
          material_mapping);

  return std::make_unique<domain::Definition<dim>>(
      std::move(mesh_ptr), finite_element_ptr);
}
template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildStamper(
    problem::ParametersI *problem_parameters,
    const std::shared_ptr<Domain> &domain_ptr,
    const std::shared_ptr<FiniteElement> &finite_element_ptr,
    const std::shared_ptr<CrossSections> &cross_sections_ptr)
-> std::unique_ptr<CFEMStamper> {

  std::unique_ptr<CFEMStamper> return_ptr = nullptr;

  // Diffusion Stamper
  if (problem_parameters->TransportModel() == problem::EquationType::kDiffusion) {

    auto diffusion_ptr = std::make_unique<formulation::scalar::CFEM_Diffusion<dim>>(
        finite_element_ptr, cross_sections_ptr);

    return_ptr = std::move(
        std::make_unique<formulation::CFEM_DiffusionStamper<dim>>(
            std::move(diffusion_ptr),
            domain_ptr,
            problem_parameters->ReflectiveBoundary()));

  } else {
    AssertThrow(false, dealii::ExcMessage("Unsuppored equation type passed"
                                         "to BuildScalarFormulation"));
  }

  return return_ptr;
}
template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildSourceUpdater(
    problem::ParametersI *,
    const std::shared_ptr<CFEMStamper> stamper_ptr)
    -> std::unique_ptr<SourceUpdater> {
  // TODO(Josh): Add option for non-gauss-seidel updating
  using SourceUpdater = iteration::updater::SourceUpdaterGaussSeidel<CFEMStamper>;
  return std::make_unique<SourceUpdater>(stamper_ptr);
}

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildParameterConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<ParameterConvergenceChecker>{

  using CheckerType = convergence::parameters::SingleParameterChecker;
  using FinalCheckerType = convergence::FinalCheckerOrN<double, CheckerType>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(single_checker_ptr));

  return_ptr->SetMaxIterations(max_iterations);

  return std::move(return_ptr);
}

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildMomentConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<MomentConvergenceChecker>{
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
auto CFEM_FrameworkBuilder<dim>::BuildFixedUpdater(
    const std::shared_ptr<CFEMStamper> &stamper_ptr)
-> std::unique_ptr<CFEM_FrameworkBuilder::FixedUpdater> {
  using FixedUpdaterType = iteration::updater::FixedUpdater<CFEMStamper>;
  return std::make_unique<FixedUpdaterType>(stamper_ptr);
}

template class CFEM_FrameworkBuilder<1>;
template class CFEM_FrameworkBuilder<2>;
template class CFEM_FrameworkBuilder<3>;



} // namespace builder

} // namespace framework

} // namespace bart