#ifndef BART_SRC_ITERATION_SUBROUTINE_TWO_GRID_HPP_
#define BART_SRC_ITERATION_SUBROUTINE_TWO_GRID_HPP_

#include "acceleration/two_grid/flux_corrector_i.hpp"
#include "calculator/residual/domain_isotropic_residual_i.hpp"
#include "framework/framework_i.hpp"
#include "iteration/subroutine/subroutine_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::iteration::subroutine {

class TwoGridAcceleration : public SubroutineI, public utility::HasDependencies {
 public:
  using FluxCorrector = acceleration::two_grid::FluxCorrectorI;
  using ResidualCalculator = calculator::residual::DomainIsotropicResidualI;
  using Framework = framework::FrameworkI;

  TwoGridAcceleration(std::unique_ptr<FluxCorrector>,
                      std::unique_ptr<Framework>,
                      std::unique_ptr<ResidualCalculator>,
                      std::shared_ptr<dealii::Vector<double>> isotropic_residual_ptr);
  auto Execute(system::System &) -> void override;

  auto flux_corrector_ptr() { return flux_corrector_ptr_.get(); };
  auto framework_ptr() { return framework_ptr_.get(); };
  auto residual_calculator_ptr() { return residual_calculator_ptr_.get(); };
  auto isotropic_residual_ptr() { return isotropic_residual_ptr_; };

 private:
  std::unique_ptr<FluxCorrector> flux_corrector_ptr_;
  std::unique_ptr<Framework> framework_ptr_;
  std::unique_ptr<ResidualCalculator> residual_calculator_ptr_;
  std::shared_ptr<dealii::Vector<double>> isotropic_residual_ptr_;
};

} // namespace bart::iteration::subroutine

#endif //BART_SRC_ITERATION_SUBROUTINE_TWO_GRID_HPP_
