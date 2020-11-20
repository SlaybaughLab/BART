#include "framework/framework_helper.hpp"

namespace bart::framework {

template<int dim>
FrameworkHelper<dim>::FrameworkHelper(const std::shared_ptr<SystemHelper>& system_helper_ptr)
    : system_helper_ptr_(system_helper_ptr) {
  AssertThrow(system_helper_ptr_ != nullptr, dealii::ExcMessage("Error in constructor of framework helper, system"
                                                                "helper pointer passed is null"))
}

template class FrameworkHelper<1>;
template class FrameworkHelper<2>;
template class FrameworkHelper<3>;

} // namespace bart::framework
