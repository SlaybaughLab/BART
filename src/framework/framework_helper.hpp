#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_

#include "framework/framework_helper_i.hpp"
#include "problem/parameters_i.h"
#include "framework/builder/framework_builder_i.hpp"
#include "material/material_protobuf.h"
#include "system/system_helper_i.hpp"
#include "system/moments/spherical_harmonic_i.h"

namespace bart::framework {

template <int dim>
class FrameworkHelper : public FrameworkHelperI<dim> {
 public:
  using SystemHelper = const system::SystemHelperI<dim>;
  FrameworkHelper(const std::shared_ptr<SystemHelper>& system_helper_ptr);

  [[nodiscard]] auto ToFrameworkParameters(const problem::ParametersI& parameters) -> framework::FrameworkParameters override;
  [[nodiscard]] auto BuildFramework(builder::FrameworkBuilderI<dim>&,
                                    framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI> override;
  [[nodiscard]] auto BuildFramework(builder::FrameworkBuilderI<dim>&,
                                    framework::FrameworkParameters&,
                                    system::moments::SphericalHarmonicI*) -> std::unique_ptr<framework::FrameworkI> override;

  auto SetSubroutineFrameworkHelper(std::unique_ptr<FrameworkHelperI<dim>> framework_helper_ptr) {
    subroutine_framework_helper_ptr_ = std::move(framework_helper_ptr); }
  auto system_helper_ptr() { return system_helper_ptr_.get(); }
  auto subroutine_framework_helper_ptr() { return subroutine_framework_helper_ptr_.get(); }
 private:
  std::shared_ptr<SystemHelper> system_helper_ptr_{ nullptr };
  // Subroune framework helper, helps build frameworks for subroutines, used only for testing
  std::unique_ptr<FrameworkHelperI<dim>> subroutine_framework_helper_ptr_{ nullptr };
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_
