#ifndef BART_DATA_SYSTEM_SYSTEM_H_
#define BART_DATA_SYSTEM_SYSTEM_H_

#include <memory>
#include <optional>

#include "system/moments/spherical_harmonic_i.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/system_types.h"
#include "system/terms/term_i.h"

namespace bart {

namespace system {

struct System {
  //TODO(Josh): Replace all calls to *_iteration_moments with *_moments

  //! Pointer to right hand side linear term
  std::unique_ptr<system::terms::MPILinearTermI> right_hand_side_ptr_;
  //! Pointer to left hand side bilinear term
  std::unique_ptr<system::terms::MPIBilinearTermI> left_hand_side_ptr_;
  //! Flux moments for the current iteration
  std::unique_ptr<system::moments::SphericalHarmonicI> current_moments = nullptr;
  //! Flux moments for the previous iteration
  std::unique_ptr<system::moments::SphericalHarmonicI> previous_moments = nullptr;
  //! System k_effective
  std::optional<double> k_effective = std::nullopt;
};

} // namespace system

} // namespace bart

#endif // BART_DATA_SYSTEM_SYSTEM_H_