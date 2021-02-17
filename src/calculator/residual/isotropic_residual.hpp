#ifndef BART_SRC_CALCULATOR_RESIDUAL_ISOTROPIC_RESIDUAL_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_ISOTROPIC_RESIDUAL_HPP_

#include <memory>

#include "calculator/residual/isotropic_residual_i.hpp"
#include "calculator/residual/vector_difference_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::calculator::residual {

class IsotropicResidual : public IsotropicResidualI, public utility::HasDependencies {
 public:
  using VectorDifferenceCalculator = calculator::residual::VectorDifferenceI;
  IsotropicResidual(std::unique_ptr<VectorDifferenceCalculator>);
  [[nodiscard]] auto CalculateIsotropicResidual(Moments *half_step_scalar_flux_ptr,
                                                Moments *previous_step_scalar_flux_ptr,
                                                const int group, const FullMatrix &sigma_s) const -> Vector override;
  auto difference_calculator_ptr() -> VectorDifferenceCalculator* { return difference_calculator_ptr_.get(); }
 private:
  std::unique_ptr<VectorDifferenceCalculator> difference_calculator_ptr_{ nullptr };
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_ISOTROPIC_RESIDUAL_HPP_
