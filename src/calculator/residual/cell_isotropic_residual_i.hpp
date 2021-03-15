#ifndef BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_I_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_I_HPP_

#include "domain/domain_types.hpp"
#include "system/moments/spherical_harmonic_i.h"

namespace bart::calculator::residual {

template <int dim>
class CellIsotropicResidualI {
 public:
  using CellPtr = typename domain::CellPtr<dim>;
  using FluxMoments = system::moments::SphericalHarmonicI;
  using Vector = dealii::Vector<double>;


  ~CellIsotropicResidualI() = default;
  virtual auto CalculateCellResidual(CellPtr, FluxMoments* current_scalar_flux_, FluxMoments* previous_scalar_flux_,
                                     int group) -> double = 0;
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_I_HPP_
