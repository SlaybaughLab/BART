#include "iteration/subroutine/two_grid_acceleration.hpp"

namespace bart::iteration::subroutine {

TwoGridAcceleration::TwoGridAcceleration(std::unique_ptr<FluxCorrector> flux_corrector_ptr,
                                         std::unique_ptr<Framework> framework_ptr,
                                         std::unique_ptr<ResidualCalculator> residual_calculator_ptr,
                                         std::shared_ptr<dealii::Vector<double>> isotropic_residual_ptr)
    : flux_corrector_ptr_(std::move(flux_corrector_ptr)),
      framework_ptr_(std::move(framework_ptr)),
      residual_calculator_ptr_(std::move(residual_calculator_ptr)),
      isotropic_residual_ptr_(std::move(isotropic_residual_ptr)) {
  std::string function_name{ "TwoGridAcceleration constructor"};
  this->AssertPointerNotNull(flux_corrector_ptr_.get(), "flux corrector", function_name);
  this->AssertPointerNotNull(framework_ptr_.get(), "framework", function_name);
  this->AssertPointerNotNull(residual_calculator_ptr_.get(), "residual calculator", function_name);
  this->AssertPointerNotNull(isotropic_residual_ptr_.get(), "isotropic residual vector", function_name);
}

auto TwoGridAcceleration::Execute(system::System& system) -> void {
  isotropic_residual_ptr_->equ(1, residual_calculator_ptr_->CalculateDomainResidual(system.current_moments.get(),
                                                                                    system.previous_moments.get()));
  framework_ptr()->SolveSystem();
  auto& error_vector = framework_ptr_->system()->current_moments->GetMoment({0, 0, 0});
  for (int group = 0; group < system.total_groups; ++group) {
    flux_corrector_ptr_->CorrectFlux(system.current_moments->GetMoment({group, 0, 0}), error_vector, group);
  }
}

} // namespace bart::iteration::subroutine