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
  std::cout << "Starting two-grid acceleration\n";
  if (!has_run_) {
    std::cout << "Running two-grid first time setup\n";
    previous_iteration_moments_ = std::make_shared<system::moments::SphericalHarmonic>(
        system.current_moments->total_groups(), 0);
    for (auto& moment : *previous_iteration_moments_) {
      moment.second.reinit(system.current_moments->begin()->second.size());
    }
    has_run_ = true;
  }

  isotropic_residual_ptr_->equ(1, residual_calculator_ptr_->CalculateDomainResidual(system.current_moments.get(),
                                                                                     previous_iteration_moments_.get()));
  std::cout << "Calculated <R>_L1 = " << isotropic_residual_ptr_->l1_norm() << '\n';

  std::cout << "Solving system\n";
  framework_ptr_->SolveSystem();
  auto error_vector = framework_ptr_->system()->current_moments->GetMoment({0, 0, 0});

  std::cout << "Correcting flux <E>_L1 = " << error_vector.l1_norm() << '\n';
  for (int group = 0; group < system.total_groups; ++group) {
    flux_corrector_ptr_->CorrectFlux(system.current_moments->GetMoment({group, 0, 0}), error_vector, group);
  }

  for (int group = 0; group < system.total_groups; ++group) {
    previous_iteration_moments_->GetMoment({group, 0, 0}) = system.current_moments->GetMoment({group, 0, 0});
  }
}

} // namespace bart::iteration::subroutine