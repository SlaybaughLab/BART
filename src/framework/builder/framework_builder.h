#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_

#include <memory>

// Problem parameters
#include "problem/parameters_i.h"

// Interface classes built by this factory
#include "quadrature/quadrature_set_i.h"
#include "convergence/reporter/mpi_i.h"

namespace bart {

namespace framework {

namespace builder {

template <int dim>
class FrameworkBuilder {
 public:
  FrameworkBuilder() = default;
  ~FrameworkBuilder() = default;

  using QuadratureSetType = quadrature::QuadratureSetI<dim>;
  using ReporterType = convergence::reporter::MpiI;

  std::shared_ptr<QuadratureSetType> BuildQuadratureSet(
      const problem::ParametersI&);

  std::unique_ptr<ReporterType> BuildConvergenceReporter();

};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
