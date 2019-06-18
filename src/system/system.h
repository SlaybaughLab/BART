#ifndef BART_DATA_SYSTEM_SYSTEM_H_
#define BART_DATA_SYSTEM_SYSTEM_H_

#include <memory>
#include <optional>

#include "system/system_types.h"
#include "data/system/term_i.h"

namespace bart {

namespace system {

struct System {
  //! Pointer to right hand side linear term
  std::unique_ptr<data::system::MPILinearTermI> right_hand_side_ptr_;
  //! Pointer to left hand side bilinear term
  std::unique_ptr<data::system::MPIBilinearTermI> left_hand_side_ptr_;
  //! Flux moments for the current iteration
  data::system::MomentsMap current_iteration_moments = {};
  //! Flux moments for the previous iteration
  data::system::MomentsMap previous_iteration_moments = {};
  //! System k_effective
  std::optional<double> k_effective = std::nullopt;
};

} // namespace system

} // namespace bart

#endif // BART_DATA_SYSTEM_SYSTEM_H_