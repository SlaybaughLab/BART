#include "calculator/residual/domain_isotropic_residual.hpp"

namespace bart::calculator::residual {

template <int dim>
DomainIsotropicResidual<dim>::DomainIsotropicResidual(
    std::unique_ptr<CellIsotropicResidualCalculator> cell_isotropic_residual_calculator_ptr,
    std::shared_ptr<Domain> domain_ptr)
    : cell_isotropic_residual_calculator_ptr_(std::move(cell_isotropic_residual_calculator_ptr)),
      domain_ptr_(std::move(domain_ptr)) {
  std::string function_name{ "DomainIsotropicResidual constructor"};
  this->AssertPointerNotNull(domain_ptr_.get(), "domain", function_name);
  this->AssertPointerNotNull(cell_isotropic_residual_calculator_ptr_.get(), "cell isotropic residual calculator",
                             function_name);
}

template<int dim>
DomainIsotropicResidualI::Vector DomainIsotropicResidual<dim>::CalculateDomainResidual(
    FluxMoments *current_flux_moments,
    FluxMoments *previous_flux_moments) {
  const int total_groups{ current_flux_moments->total_groups() };
  AssertThrow(total_groups == previous_flux_moments->total_groups(),
              dealii::ExcMessage("Error in CalculateDomainResidual, flux moments total groups inequal"))
  Vector isotropic_residual(this->domain_ptr_->total_degrees_of_freedom());
  auto cell_vector = this->domain_ptr_->GetCellVector();
  std::vector<unsigned int> cell_global_dofs_indices(cell_vector.size());

  for (auto& cell : domain_ptr_->Cells()) {
    cell_vector = 0;
    cell->get_dof_indices(cell_global_dofs_indices);
    for (int group = 0; group < total_groups; ++group) {
      cell_vector.add(this->cell_isotropic_residual_calculator_ptr_->CalculateCellResidual(cell, current_flux_moments,
                                                                                           previous_flux_moments,
                                                                                           group));
    }
    isotropic_residual.add(cell_global_dofs_indices, cell_vector);
  }

  return isotropic_residual;
}

template class DomainIsotropicResidual<1>;
template class DomainIsotropicResidual<2>;
template class DomainIsotropicResidual<3>;

} // namespace bart::calculator::residual
