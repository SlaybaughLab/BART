#ifndef BART_SRC_ITERATION_UPDATER_FIXED_UDPATER_I_H_
#define BART_SRC_ITERATION_UPDATER_FIXED_UDPATER_I_H_

#include "system/system.h"
#include "system/system_types.h"

namespace bart {

namespace data {

struct System;

} // namespace data

namespace iteration {

namespace updater {
/*! \brief Interface for classes that update fixed terms in systems.
 */
class [[deprecated]] FixedUpdaterI {
 public:
  virtual ~FixedUpdaterI() = default;
  /*! \brief Updates the fixed term in the provided system.
   *
   * Depending on implementation, this function will access the bilinear (left
   * hand side), linear (right hand side), or both terms and update them with
   * the appropriate fixed terms. The terms may depend on group, angle, or both.
   *
   * @param system system to update.
   * @param group group to update.
   * @param angle angle to update.
   */
  virtual void UpdateFixedTerms(system::System& system,
                                system::GroupNumber group,
                                system::AngleIndex angle) = 0;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_FIXED_UDPATER_I_H_