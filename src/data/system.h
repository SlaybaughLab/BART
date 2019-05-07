#ifndef BART_DATA_SYSTEM_SYSTEM_H_
#define BART_DATA_SYSTEM_SYSTEM_H_

#include <memory>

#include "data/system/system_types.h"
#include "data/system/right_hand_side_i.h"

namespace bart {

namespace data {

struct System {
  std::unique_ptr<system::RightHandSideI> right_hand_side_ptr_;
};
} // namespace data

} // namespace bart

#endif // BART_DATA_SYSTEM_SYSTEM_H_