#ifndef BART_SRC_ACCELERATION_TWO_GRID_FLUX_CORRECTOR_I_HPP_
#define BART_SRC_ACCELERATION_TWO_GRID_FLUX_CORRECTOR_I_HPP_

#include <deal.II/lac/vector.h>

#include "utility/has_description.h"

namespace bart::acceleration::two_grid {

class FluxCorrectorI : public utility::HasDescription{
 public:
  using Vector = dealii::Vector<double>;
  virtual ~FluxCorrectorI() = default;
  virtual auto CorrectFlux(Vector& flux_to_correct, const Vector& error, const int group) const -> void = 0;
};

} // namespace bart::acceleration::two_grid

#endif //BART_SRC_ACCELERATION_TWO_GRID_FLUX_CORRECTOR_I_HPP_
