#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_

#include "problem/parameters_i.h"
#include "framework/builder/framework_builder_i.hpp"
#include "material/material_protobuf.h"
#include "system/system_helper_i.hpp"
#include "system/moments/spherical_harmonic_i.h"

namespace bart::framework {

template <int dim>
class FrameworkHelper {
 public:
  using SystemHelper = const system::SystemHelperI<dim>;
  FrameworkHelper(const std::shared_ptr<SystemHelper>& system_helper_ptr);

  [[nodiscard]] auto ToFrameworkParameters(const problem::ParametersI& parameters) -> framework::FrameworkParameters;
  [[nodiscard]] auto BuildFramework(builder::FrameworkBuilderI<dim>&,
                                    framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI>;
  [[nodiscard]] auto BuildFramework(builder::FrameworkBuilderI<dim>&,
                                    framework::FrameworkParameters&,
                                    system::moments::SphericalHarmonicI*) -> std::unique_ptr<framework::FrameworkI>;

  auto system_helper_ptr() { return system_helper_ptr_.get(); }
 private:
  std::shared_ptr<SystemHelper> system_helper_ptr_{ nullptr };
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_
