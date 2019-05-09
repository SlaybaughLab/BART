#include "iteration/updater/source_updater.h"
#include "formulation/cfem_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
std::shared_ptr<data::system::MPIVector> SourceUpdater<StamperType>::GetSourceVectorPtr(
    VariableTerms term,
    data::System& system,
    data::system::GroupNumber group,
    data::system::AngleIndex angle) {
  auto source_vector_ptr =
      system.right_hand_side_ptr_->GetVariableTermPtr({group, angle}, term);

  if (source_vector_ptr == nullptr) {
    std::ostringstream oss;
    oss << "Right hand sidse returned nullptr for group " << group << " angle"
        << angle << "combination";
    AssertThrow(false, dealii::ExcMessage(oss.str()));
  }
  return source_vector_ptr;
}

template class SourceUpdater<formulation::CFEMStamperI>;

} // namespace updater

} // namespace iteration

} // namespace bart