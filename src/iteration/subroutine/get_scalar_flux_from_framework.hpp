#ifndef BART_SRC_ITERATION_SUBROUTINE_GET_SCALAR_FLUX_FROM_FRAMEWORK_HPP_
#define BART_SRC_ITERATION_SUBROUTINE_GET_SCALAR_FLUX_FROM_FRAMEWORK_HPP_

#include "framework/framework_i.hpp"
#include "iteration/subroutine/subroutine_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::iteration::subroutine {

class GetScalarFluxFromFramework : public SubroutineI, public utility::HasDependencies {
 public:
  using Framework = bart::framework::FrameworkI;
  GetScalarFluxFromFramework(std::unique_ptr<Framework> framework_ptr)
      : framework_ptr_(std::move(framework_ptr)) {
    AssertPointerNotNull(framework_ptr_.get(), "framework pointer", "GetScalarFluxFromFramework constructor");
  }

  auto Execute(system::System& system) -> void override {
    this->framework_ptr_->SolveSystem();
    auto& solved_scalar_fluxes = *this->framework_ptr_->system()->current_moments;
    for (int group = 0; group < system.total_groups; ++group) {
      std::array<int, 3> index{group, 0, 0};
      auto& solved_scalar_flux = solved_scalar_fluxes.GetMoment(index);
      auto& system_scalar_flux = system.current_moments->GetMoment(index);
      system_scalar_flux = solved_scalar_flux;
    }
  };

  auto framework_ptr() -> Framework* { return framework_ptr_.get(); }
 private:
  std::unique_ptr<Framework> framework_ptr_{ nullptr };
};

} // namespace bart::iteration::subroutine

#endif //BART_SRC_ITERATION_SUBROUTINE_GET_SCALAR_FLUX_FROM_FRAMEWORK_HPP_
