#include "formulation/updater/saaf_updater.h"

namespace bart {

namespace formulation {

namespace updater {

template<int dim>
SAAFUpdater<dim>::SAAFUpdater(std::unique_ptr<SAAFFormulationType> formulation_ptr,
                              std::unique_ptr<StamperType> stamper_ptr)
    : formulation_ptr_(std::move(formulation_ptr)),
      stamper_ptr_(std::move(stamper_ptr)) {
  AssertThrow(formulation_ptr_ != nullptr,
      dealii::ExcMessage("Error in constructor of SAAFUpdater, formulation "
                         "pointer passed is null"))
  AssertThrow(stamper_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of SAAFUpdater, stamper "
                                 "pointer passed is null"))
}

template class SAAFUpdater<1>;
template class SAAFUpdater<2>;
template class SAAFUpdater<3>;

} // namespace updater

} // namespace formulation

} // namespace bart
