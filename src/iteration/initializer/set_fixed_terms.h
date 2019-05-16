#ifndef BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_
#define BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_

#include <memory>

#include "data/system.h"
#include "iteration/initializer/initializer_i.h"

namespace bart {

//namespace data { struct System; } // namespace data

namespace iteration {

namespace updater { class FixedUpdaterI;} // namespace updater

namespace initializer {
/*! \brief Initializes a system by setting the fixed terms.
 *
 * This class takes a class of type iteration::updater::FixedUpdaterI and
 * uses the UpdateFixedTerms method to initialize the system.
 *
 */
class SetFixedTerms : public InitializerI {
 public:
  /*! \brief Constructor, takes ownership of fixed updater.
   *
   * Move semantics must be used in construction.
   *
   * @param fixed_updater_ptr fixed updater that will be owned by this
   *                          initializer.
   * @param total_groups total energy groups
   * @param total_angles total angles
   */
  SetFixedTerms(std::unique_ptr<updater::FixedUpdaterI> fixed_updater_ptr,
                const int total_groups,
                const int total_angles);

  virtual ~SetFixedTerms() = default;
  virtual void Initialize(data::System& sys) override;

  /*! \brief Get pointer to the owned fixed updater.
   *
   * @return pointer to owned fixed updater.
   */
  updater::FixedUpdaterI* fixed_updater_ptr() const {
    return fixed_updater_ptr_.get(); };

  /*! \brief Get total energy groups.
   *
   * @return total energy groups.
   */
  int total_groups() const { return total_groups_; }
  /*! \brief Get total angles.
   *
   * @return total angles.
   */
  int total_angles() const { return total_angles_; }

 protected:
  std::unique_ptr<updater::FixedUpdaterI> fixed_updater_ptr_;
  const int total_groups_;
  const int total_angles_;
};

} // namespace initializer

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_INITIALIZER_SET_FIXED_TERMS_H_