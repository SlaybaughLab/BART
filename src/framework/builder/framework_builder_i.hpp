#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_

#include <memory>
#include <string>

#include "framework/framework_i.hpp"
#include "framework/framework_parameters.hpp"

namespace bart::framework::builder {

template <int dim>
class FrameworkBuilderI {
 public:
  using FrameworkI = framework::FrameworkI;
};

template <int dim>
auto BuildFramework(FrameworkBuilderI<dim>&,
                    const framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI>;

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_
