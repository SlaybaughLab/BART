#include "iteration/updater/fixed_updater.h"
#include "formulation/cfem_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
FixedUpdater<StamperType>::FixedUpdater(
    std::unique_ptr<StamperType> stamper_ptr) {
  AssertThrow(stamper_ptr != nullptr,
      dealii::ExcMessage("Error in constructor of FixedUpdater: "
                         "stamper pointer is null."));
  stamper_ptr_ = std::move(stamper_ptr);
}

template <typename StamperType>
void FixedUpdater<StamperType>::UpdateFixedTerms(data::System& system,
                                                 data::system::GroupNumber group,
                                                 data::system::AngleIndex angle) {

}

template class FixedUpdater<formulation::CFEMStamperI>;

} // namespace updater

} // namespace iteration

} // namespace bart