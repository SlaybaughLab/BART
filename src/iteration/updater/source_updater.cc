#include <formulation/angular_stamper_i.h>
#include "iteration/updater/source_updater.h"
#include "formulation/cfem_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
std::shared_ptr<system::MPIVector> SourceUpdater<StamperType>::GetSourceVectorPtr(
    VariableTerms term,
    system::System& system,
    system::GroupNumber group,
    system::AngleIndex angle) {
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
template class SourceUpdater<formulation::AngularStamperI<1>>;
template class SourceUpdater<formulation::AngularStamperI<2>>;
template class SourceUpdater<formulation::AngularStamperI<3>>;

} // namespace updater

} // namespace iteration

} // namespace bart