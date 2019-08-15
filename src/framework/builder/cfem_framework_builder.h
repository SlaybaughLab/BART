#ifndef BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_

#include <memory>

#include "data/cross_sections.h"
#include "domain/definition_i.h"
#include "domain/finite_element/finite_element_i.h"
#include "problem/parameters_i.h"
#include "framework/builder/framework_builder_i.h"
#include "formulation/cfem_stamper_i.h"
#include "iteration/updater/source_updater_i.h"
#include "iteration/updater/fixed_updater_i.h"
#include "convergence/final_i.h"
#include "solver/group/single_group_solver_i.h"
#include "system/system_types.h"


namespace bart {

namespace framework {

namespace builder {

template <int dim>
class CFEM_FrameworkBuilder : public FrameworkBuilderI {
 public:
  using CrossSections = data::CrossSections;
  using Domain = typename domain::DefinitionI<dim>;
  using CFEMStamper = formulation::CFEMStamperI;
  using FiniteElement = typename domain::finite_element::FiniteElementI<dim>;
  using FixedUpdater = iteration::updater::FixedUpdaterI;
  using SingleGroupSolver = solver::group::SingleGroupSolverI;
  using SourceUpdater = iteration::updater::SourceUpdaterI;
  using ParameterConvergenceChecker = convergence::FinalI<double>;
  using MomentConvergenceChecker = convergence::FinalI<system::moments::MomentVector>;

  CFEM_FrameworkBuilder() = default;
  virtual ~CFEM_FrameworkBuilder() = default;

  /*! \brief Returns a FiniteElement object.
   *
   * Only FiniteElementGaussian is implemented. Future versions may include a
   * problem parameter that specifies the TYPE of finite element object.
   *
   * @param problem_parameters problem parameters
   * @return
   */
  std::unique_ptr<FiniteElement> BuildFiniteElement(
      problem::ParametersI* problem_parameters);

  std::unique_ptr<Domain> BuildDomain(
      problem::ParametersI* problem_parameters,
      const std::shared_ptr<FiniteElement> &finite_element_ptr,
      std::string material_mapping);

  std::unique_ptr<FixedUpdater> BuildFixedUpdater(
      const std::shared_ptr<CFEMStamper> &stamper_ptr);

  std::unique_ptr<SingleGroupSolver> BuildSingleGroupSolver(
      const int max_iterations = 100,
      const double convergence_tolerance = 1e-10);

  std::unique_ptr<CFEMStamper> BuildStamper(
      problem::ParametersI* problem_parameters,
      const std::shared_ptr<Domain> &domain_ptr,
      const std::shared_ptr<FiniteElement> &finite_element_ptr,
      const std::shared_ptr<CrossSections> &cross_sections_ptr);

  std::unique_ptr<SourceUpdater> BuildSourceUpdater(
      problem::ParametersI* problem_parameters,
      const std::shared_ptr<CFEMStamper> stamper_ptr);

  std::unique_ptr<ParameterConvergenceChecker> BuildParameterConvergenceChecker(
      double max_delta, int max_iterations);

  std::unique_ptr<MomentConvergenceChecker> BuildMomentConvergenceChecker(
      double max_delta, int max_iterations);

};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_
