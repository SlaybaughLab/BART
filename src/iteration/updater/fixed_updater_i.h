#ifndef BART_SRC_ITERATION_UPDATER_FIXED_UDPATER_I_H_
#define BART_SRC_ITERATION_UPDATER_FIXED_UDPATER_I_H_

#include "data/system/system_types.h"

namespace bart {

namespace data {

struct System;

} // namespace data

namespace iteration {

namespace updater {

class FixedUpdaterI {
 public:
  virtual ~FixedUpdaterI() = default;
  virtual void UpdateFixedTerms(data::System& system,
                                data::system::GroupNumber group,
                                data::system::AngleIndex angle) = 0;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_FIXED_UDPATER_I_H_