#ifndef BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_

#include <memory>

#include "problem/parameters_i.h"
#include "formulation/cfem_stamper_i.h"
#include "framework/builder/framework_builder_i.h"


namespace bart {

namespace framework {

namespace builder {

template <int dim>
class CFEM_FrameworkBuilder : public FrameworkBuilderI {
 public:
  CFEM_FrameworkBuilder() = default;
  virtual ~CFEM_FrameworkBuilder() = default;

  std::shared_ptr<formulation::CFEMStamperI> BuildStamper(
      problem::ParametersI* problem_parameters,
      std::string material_mapping);
};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_
