#include "angular_source_updater_gauss_seidel.h"

#include "formulation/angular_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template<typename StamperType>
AngularSourceUpdaterGaussSeidel<StamperType>::AngularSourceUpdaterGaussSeidel(
    std::shared_ptr<StamperType> stamper_ptr,
    std::shared_ptr<QuadratureSetType> quadrature_set_ptr)
    : SourceUpdater<StamperType>(stamper_ptr),
      quadrature_set_ptr_(quadrature_set_ptr) {
  AssertThrow(this->stamper_ptr_ != nullptr,
      dealii::ExcMessage("Error in constructor of AngularSourceUpdaterGaussSeidel,"
                         " stamper pointer is null"))
  AssertThrow(quadrature_set_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of AngularSourceUpdaterGaussSeidel,"
                                 " quadrature set pointer is null"))
}

template<typename StamperType>
void AngularSourceUpdaterGaussSeidel<StamperType>::UpdateScatteringSource(
    system::System &system,
    system::GroupNumber group,
    system::AngleIndex angle) {
  using term = system::terms::VariableLinearTerms;
  auto scattering_source_vector_ptr_ = this->GetSourceVectorPtr(
      term::kScatteringSource, system, group, angle);

  *scattering_source_vector_ptr_ = 0;

  const auto& moments = system.current_moments->moments();
  const auto quadrature_point = this->quadrature_set_ptr_->GetQuadraturePoint(
      quadrature::QuadraturePointIndex(angle));

  this->stamper_ptr_->StampScatteringSourceTerm(*scattering_source_vector_ptr_,
                                                quadrature_point,
                                                system::EnergyGroup(group),
                                                moments.at({group, 0, 0}),
                                                moments);
}
template<typename StamperType>
void AngularSourceUpdaterGaussSeidel<StamperType>::UpdateFissionSource(
    system::System &system,
    system::GroupNumber group,
    system::AngleIndex angle) {
  const std::optional<double> k_eff = system.k_effective;
  AssertThrow(k_eff.has_value(),
              dealii::ExcMessage("Error in UpdateFissionSource, system.k_effective"
                                 " is invalid (has no value)."))
  AssertThrow(k_eff.value() > 0,
              dealii::ExcMessage("Error in UpdateFissionSource, system.k_effective"
                                 " is invalid (less than or equal to 0)."))

  using term = system::terms::VariableLinearTerms;
  auto fission_source_vector_ptr_ = this->GetSourceVectorPtr(
      term::kFissionSource, system, group, angle);

  *fission_source_vector_ptr_ = 0;

  const auto& moments = system.current_moments->moments();
  const auto quadrature_point = this->quadrature_set_ptr_->GetQuadraturePoint(
      quadrature::QuadraturePointIndex(angle));

  this->stamper_ptr_->StampFissionSourceTerm(*fission_source_vector_ptr_,
                                             quadrature_point,
                                             system::EnergyGroup(group),
                                             system.k_effective.value(),
                                             moments.at({group, 0, 0}),
                                             moments);
}

template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<1>>;
template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<2>>;
template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<3>>;

} // namespace updater

} // namespace iteration

} // namespace bart
