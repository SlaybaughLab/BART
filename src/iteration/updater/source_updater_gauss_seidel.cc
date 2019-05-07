#include "iteration/updater/source_updater_gauss_seidel.h"

namespace bart {

namespace iteration {

namespace updater {

template <>
void SourceUpdaterGaussSeidel<formulation::CFEMStamperI>::UpdateScatteringSource(
    data::System& system,
    data::system::GroupNumber group,
    data::system::AngleIndex angle) {

  using VariableTerms = data::system::RightHandSideI::VariableTerms;

  auto scattering_source_vector_ptr_ =
      system.right_hand_side_ptr_->GetVariablePtr({group, angle},
                                                  VariableTerms::kScatteringSource);

  data::MomentVector& in_group_moment = system.current_iteration_moments[{group, 0, 0}];
  data::MomentsMap& out_group_moments = system.current_iteration_moments;

  *scattering_source_vector_ptr_ = 0;

  stamper_ptr_->StampScatteringSource(*scattering_source_vector_ptr_,
                                      group,
                                      in_group_moment,
                                      out_group_moments);
}

template class SourceUpdaterGaussSeidel<formulation::CFEMStamperI>;

} // namespace updater

} // namespace iteration

} // namespace bart