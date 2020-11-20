#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_

#include "framework/builder/framework_builder_i.hpp"
#include "system/system_helper_i.hpp"

namespace bart::framework {

template <int dim>
class FrameworkHelper {
 public:
  using SystemHelper = const system::SystemHelperI<dim>;
  FrameworkHelper(const std::shared_ptr<SystemHelper>& system_helper_ptr);

  [[nodiscard]] auto BuildFramework(builder::FrameworkBuilderI<dim>&,
                                    const framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI>;

  auto system_helper_ptr() { return system_helper_ptr_.get(); }
 private:
  std::shared_ptr<SystemHelper> system_helper_ptr_{ nullptr };
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_HELPER_HPP_
