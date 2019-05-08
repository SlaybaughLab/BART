#ifndef BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_
#define BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_

#include <memory>

#include "data/system/system_types.h"
#include "iteration/updater/source_updater_i.h"
#include "formulation/cfem_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
class SourceUpdater : public SourceUpdaterI {
 public:
  using VariableTerms = data::system::RightHandSideI::VariableTerms;
  using MPIVector = data::system::MPIVector;

  explicit SourceUpdater(std::unique_ptr<StamperType> stamper_ptr)
      : stamper_ptr_(std::move(stamper_ptr)) {};
  virtual ~SourceUpdater() override = default;

  std::shared_ptr<MPIVector> GetSourceVectorPtr(data::System& system,
                                                data::system::GroupNumber group,
                                                data::system::AngleIndex angle,
                                                VariableTerms term) {
    auto source_vector_ptr =
        system.right_hand_side_ptr_->GetVariablePtr({group, angle}, term);

    if (source_vector_ptr == nullptr) {
      std::ostringstream oss;
      oss << "Right hand sidse returned nullptr for group " << group << " angle"
          << angle << "combination";
      AssertThrow(false, dealii::ExcMessage(oss.str()));
    }
    return source_vector_ptr;
  }

 protected:
  std::unique_ptr<StamperType> stamper_ptr_;

};

} // namespace updater

} // namespace iteration

} // namespace bart

#endif // BART_SRC_ITERATION_UPDATER_SOURCE_UPDATER_H_