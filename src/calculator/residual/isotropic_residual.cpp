#include "calculator/residual/isotropic_residual.hpp"

namespace bart::calculator::residual {

IsotropicResidual::IsotropicResidual(std::unique_ptr<VectorDifferenceCalculator> vector_difference_calculator_ptr)
: difference_calculator_ptr_(std::move(vector_difference_calculator_ptr)) {
  AssertPointerNotNull(difference_calculator_ptr_.get(), "vector difference calculator", "IsotropicResidual const'r");
}

auto IsotropicResidual::CalculateIsotropicResidual(Moments *half_step_scalar_flux_ptr,
                                                   Moments *previous_step_scalar_flux_ptr,
                                                   const int group,
                                                   const FullMatrix &sigma_s) const -> Vector {
  // Quality check inputs
  const std::string error_text{ "Error in IsotropicResidual::CalculateIsotropicResidual: "};
  AssertThrow(group >= 0, dealii::ExcMessage(error_text + "group must be a positive value"))
  const int total_groups = half_step_scalar_flux_ptr->total_groups();
  AssertThrow(total_groups == previous_step_scalar_flux_ptr->total_groups(),
              dealii::ExcMessage(error_text + "half step groups and previous step groups mismatched."))
  AssertThrow(group < total_groups, dealii::ExcMessage(error_text + "specified group must be < total_groups"))


  Vector return_vector;

  for (int group_in = group + 1; group_in < total_groups; ++group_in) {
    std::array<int, 3> index{ group_in, 0, 0};
    auto& half_step_scalar_flux = half_step_scalar_flux_ptr->GetMoment(index);
    auto& previous_step_scalar_flux = previous_step_scalar_flux_ptr->GetMoment(index);
    auto group_residual = difference_calculator_ptr_->CalculateResidual(half_step_scalar_flux,
                                                                        previous_step_scalar_flux,
                                                                        sigma_s(group, group_in));
    if (return_vector.size() == 0) {
      return_vector = group_residual;
    } else {
      return_vector.add(1, group_residual);
    }
  }
  return return_vector;
}

} // namespace bart::calculator::residual
