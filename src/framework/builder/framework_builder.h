#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_

#include <memory>

// Problem parameters
#include "problem/parameters_i.h"

// Interface classes built by this factory
#include "convergence/reporter/mpi_i.h"
#include "domain/finite_element/finite_element_i.h"
#include "quadrature/quadrature_set_i.h"


namespace bart {

namespace framework {

namespace builder {

template <int dim>
class FrameworkBuilder {
 public:
  FrameworkBuilder() = default;
  ~FrameworkBuilder() = default;

  using ParametersType = const problem::ParametersI&;

  using FiniteElementType = domain::finite_element::FiniteElementI<dim>;
  using QuadratureSetType = quadrature::QuadratureSetI<dim>;
  using ReporterType = convergence::reporter::MpiI;

  std::unique_ptr<ReporterType> BuildConvergenceReporter();
  std::unique_ptr<FiniteElementType> BuildFiniteElement(ParametersType);
  std::shared_ptr<QuadratureSetType> BuildQuadratureSet(ParametersType);
};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
