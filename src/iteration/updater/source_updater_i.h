#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_

#include "system/system_types.h"
#include "system/system.h"

namespace bart {

namespace iteration {

namespace updater {

/*! \brief Interface for classes that update source terms in systems.
 */
class SourceUpdaterI {
 public:
  virtual ~SourceUpdaterI() = default;
   /*! \brief Updates the scattering source term in the provided system.
   *
   * Depending on implementation, this function will access the bilinear (left
   * hand side), linear (right hand side), or both terms and update them with
   * the appropriate fixed terms. The terms may depend on group, angle, or both.
   *
   * @param system system to update.
   * @param group group to update.
   * @param angle angle to update.
   */
  virtual void UpdateScatteringSource(system::System& system,
                                      data::system::GroupNumber group,
                                      data::system::AngleIndex angle) = 0;
  /*! \brief Updates the fission source term in the provided system.
   *
   * Depending on implementation, this function will access the bilinear (left
   * hand side), linear (right hand side), or both terms and update them with
   * the appropriate fixed terms. The terms may depend on group, angle, or both.
   *
   * @param system system to update.
   * @param group group to update.
   * @param angle angle to update.
   */
  virtual void UpdateFissionSource(system::System& system,
                                   data::system::GroupNumber group,
                                   data::system::AngleIndex angle) = 0;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_I_H_