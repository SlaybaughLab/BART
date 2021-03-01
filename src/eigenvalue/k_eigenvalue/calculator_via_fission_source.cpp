#include "eigenvalue/k_eigenvalue/calculator_via_fission_source.hpp"

namespace bart::eigenvalue::k_eigenvalue {

CalculatorViaFissionSource::CalculatorViaFissionSource(
    std::unique_ptr<FissionSourceCalculator> fission_source_calculator_ptr, const double initial_k_eigenvalue,
    const double initial_fission_source)
    : fission_source_calculator_(std::move(fission_source_calculator_ptr)),
      initial_k_eigenvalue_(initial_k_eigenvalue),
      initial_fission_source_(initial_fission_source) {
  AssertThrow(initial_k_eigenvalue > 0.0, dealii::ExcMessage("Error in constructor of CalculatorViaFissionSource, "
                                                             "initial k_eigenvalue must be > 0"));
  AssertThrow(initial_fission_source > 0.0, dealii::ExcMessage("Error in constructor of CalculatorViaFissionSource, "
                                                               "initial fission source must be > 0"));
  this->set_description("k-eigenvalue updater via fission source.");
}

double CalculatorViaFissionSource::CalculateK_Eigenvalue(system::System &system) {
  current_fission_source_ = fission_source_calculator_->AggregatedFissionSource(system.current_moments.get());
  AssertThrow(current_fission_source_ > 0, dealii::ExcMessage("Error in CalculateK_Eigenvalue, fission source is 0"));
  k_eigenvalue_ = initial_k_eigenvalue_ * current_fission_source_.value() / initial_fission_source_;
  return k_eigenvalue_.value();
}

} // namespace bart::eigenvalue::k_eigenvalue