#include "formulation/updater/diffusion_updater.h"

namespace bart {

namespace formulation {

namespace updater {

template<int dim>
DiffusionUpdater<dim>::DiffusionUpdater(
    std::unique_ptr<DiffusionFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr)
    : formulation_ptr_(std::move(formulation_ptr)),
      stamper_ptr_(std::move(stamper_ptr)) {
  AssertThrow(formulation_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of DiffusionUpdater, "
                                 "formulation pointer passed is null"))
  AssertThrow(stamper_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of DiffusionUpdater, "
                                 "stamper pointer passed is null"))
}

template<int dim>
DiffusionUpdater<dim>::DiffusionUpdater(
    std::unique_ptr<DiffusionFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr,
    std::unordered_set<problem::Boundary> reflective_boundaries)
    : DiffusionUpdater(std::move(formulation_ptr), std::move(stamper_ptr)) {
  reflective_boundaries_ = reflective_boundaries;
}

template class DiffusionUpdater<1>;
template class DiffusionUpdater<2>;
template class DiffusionUpdater<3>;

} // namespace updater

} // namespace formulation

} // namespace bart
