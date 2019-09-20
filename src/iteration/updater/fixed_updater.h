#ifndef BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_
#define BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_

#include <memory>

#include "iteration/updater/fixed_updater_i.h"

namespace bart {

namespace iteration {

namespace updater {

/*! \brief Updates the fixed term in a system using a provided stamper.
 *
 * Implementation of this class must be created for each stamper, as they may
 * have different fixed terms. The entire system is passed to to the UpdateFixedTerms
 * method, allowing for updating of any part of the system.
 *
 * @tparam StamperType the type of stamper that this class will use to update
 *                     the system.
 */
template <typename StamperType>
class FixedUpdater : public FixedUpdaterI {
 public:
  /*! \brief Constructor, takes ownership of stamper.
   *
   * The stamper is passed via a shared pointer, so when constructing, move
   * semantics must be used.
   *
   * Example:
   * \code
   * auto stamper_ptr = std::make_shared<MyStamperType>();
   * FixedUpdater<MyStamperType> updater(std::move(stamper_ptr)); // note use of std::move
   * \endcode
   *
   * @param stamper_ptr shared pointer to the stamper that this class will take
   *                    ownership of.
   */
  explicit FixedUpdater(std::shared_ptr<StamperType> stamper_ptr);
  /*! \brief Destructor.
   *
   * Marked virtual to allow deriving from this class.
   *
   */
  virtual ~FixedUpdater() = default;
  void UpdateFixedTerms(system::System& system,
                        system::GroupNumber group,
                        system::AngleIndex angle) override;
  /*! \brief Returns a pointer to the stored stamper.
   *
   * @return a raw pointer to the stamper.
   */
  StamperType* GetStamperPtr() const {
    return stamper_ptr_.get();
  };

 protected:
  //! Stored pointer to dependent stamper
  std::shared_ptr<StamperType> stamper_ptr_;
};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_FIXED_UPDATER_H_