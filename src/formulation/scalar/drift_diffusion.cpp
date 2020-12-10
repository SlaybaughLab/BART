#include "formulation/scalar/drift_diffusion.hpp"

namespace bart::formulation::scalar {

template<int dim>
DriftDiffusion<dim>::DriftDiffusion(std::shared_ptr<FiniteElement> finite_element_ptr,
                                    std::shared_ptr<CrossSections> cross_sections_ptr,
                                    std::shared_ptr<DriftDiffusionCalculator> drift_diffusion_calculator_ptr)
    : finite_element_ptr_(finite_element_ptr),
      cross_sections_ptr_(cross_sections_ptr),
      drift_diffusion_calculator_ptr_(drift_diffusion_calculator_ptr) {
  std::string function_name{"formulation::scalar::DriftDiffusion constructor"};
  AssertPointerNotNull(finite_element_ptr_.get(), "finite_element_ptr", function_name);
  AssertPointerNotNull(cross_sections_ptr_.get(), "cross_sections_ptr", function_name);
  AssertPointerNotNull(drift_diffusion_calculator_ptr_.get(), "drift_diffusion_calculator_ptr", function_name);
  cell_quadrature_points_ = finite_element_ptr->n_cell_quad_pts();
  cell_degrees_of_freedom_ = finite_element_ptr->dofs_per_cell();
}
template<int dim>
auto DriftDiffusion<dim>::FillCellDriftDiffusionTerm(Matrix& to_fill,
                                                     const CellPtr& cell_ptr,
                                                     system::EnergyGroup group,
                                                     const Vector& group_scalar_flux,
                                                     const Vector& integrated_angular_flux) const -> void {
  finite_element_ptr_->SetCell(cell_ptr);
  const int material_id{ cell_ptr->material_id() };
  const double sigma_t{ cross_sections_ptr_->sigma_t.at(material_id).at(group.get()) };
  const double diffusion_coeff{ cross_sections_ptr_->diffusion_coef.at(material_id).at(group.get()) };

  const auto scalar_flux_at_q{ finite_element_ptr_->ValueAtQuadrature(group_scalar_flux) };
  const auto integrated_angular_flux_at_q{ finite_element_ptr_->ValueAtQuadrature(integrated_angular_flux) };

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian{ finite_element_ptr_->Jacobian(q) };
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      const auto drift_diffusion{drift_diffusion_calculator_ptr_->DriftDiffusion(
          scalar_flux_at_q.at(q),
          integrated_angular_flux_at_q.at(q),
          finite_element_ptr_->ShapeGradient(i, q),
          sigma_t,
          diffusion_coeff)};
      const auto drift_diffusion_term{ drift_diffusion * finite_element_ptr_->ShapeValue(i, q) };
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i,j) += drift_diffusion_term * this->finite_element_ptr_->ShapeGradient(j, q) * jacobian;
      }
    }
  }
}

template class DriftDiffusion<1>;
template class DriftDiffusion<2>;
template class DriftDiffusion<3>;

} // namespace bart::formulation::scalar
