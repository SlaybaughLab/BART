#ifndef BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_

#include <memory>

#include "domain/finite_element_i.h"
#include "problem/parameters_i.h"
#include "formulation/cfem_stamper_i.h"
#include "framework/builder/framework_builder_i.h"


namespace bart {

namespace framework {

namespace builder {

template <int dim>
class CFEM_FrameworkBuilder : public FrameworkBuilderI {
 public:
  using FiniteElement = typename domain::FiniteElementI<dim>;

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
  std::shared_ptr<FiniteElement> BuildFiniteElement(
      problem::ParametersI* problem_parameters);

  std::shared_ptr<formulation::CFEMStamperI> BuildStamper(
      problem::ParametersI* problem_parameters,
      std::string material_mapping);
};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_CFEM_FRAMEWORK_BUILDER_H_
