#include "framework/builder/framework_builder_i.hpp"

namespace bart::framework::builder {

template <int dim>
auto BuildFramework(FrameworkBuilderI<dim>&,
                    const framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI> {
  return nullptr;
}

template auto BuildFramework(FrameworkBuilderI<1>&, const framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI>;
template auto BuildFramework(FrameworkBuilderI<2>&, const framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI>;
template auto BuildFramework(FrameworkBuilderI<3>&, const framework::FrameworkParameters&) -> std::unique_ptr<framework::FrameworkI>;

} // namespace bart::framework::builder
