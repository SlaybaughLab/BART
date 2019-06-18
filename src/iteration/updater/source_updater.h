#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_

#include <memory>

#include "system/system_types.h"
#include "system/terms/term_types.h"
#include "iteration/updater/source_updater_i.h"
#include "formulation/cfem_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

/*! \brief Updates the source terms in a system using a provided stamper.
 *
 * Implementation of this class must be created for each stamper, as they may
 * have different methods for updating source terms. The entire system is passed
 * to the member functions, allowing for updating or accessing of any part of
 * the system.
 *
 * @tparam StamperType the type of stamper that this class will use to update
 *                     the system.
 */
template <typename StamperType>
class SourceUpdater : public SourceUpdaterI {
 public:
  using VariableTerms = system::terms::VariableLinearTerms;
  using MPIVector = system::MPIVector;

  /*! \brief Constructor, takes ownership of stamper.
   *
   * The stamper is passed via a unique pointer, so when constructing, move
   * semantics must be used.
   *
   * Example:
   * \code
   * auto stamper_ptr = std::make_unique<MyStamperType>();
   * SourceUpdater<MyStamperType> updater(std::move(stamper_ptr)); // note use of std::move
   * \endcode
   *
   * @param stamper_ptr unique pointer to the stamper that this class will take
   *                    ownership of.
   */
  explicit SourceUpdater(std::unique_ptr<StamperType> stamper_ptr)
      : stamper_ptr_(std::move(stamper_ptr)) {};
  /*! \brief Destructor.
   *
   * Marked virtual to allow deriving from this class.
   *
   */
  virtual ~SourceUpdater() override = default;

  /*! Returns the right hand side vector from the system for a specific source term.
   *
   * @param term source term
   * @param system system holding the right hand side data
   * @param group right hand side group
   * @param angle right hand side angle
   * @return pointer to the right hand side
   */
  std::shared_ptr<MPIVector> GetSourceVectorPtr(VariableTerms term,
                                                system::System& system,
                                                system::GroupNumber group,
                                                system::AngleIndex angle);
 protected:
  std::unique_ptr<StamperType> stamper_ptr_;

};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_