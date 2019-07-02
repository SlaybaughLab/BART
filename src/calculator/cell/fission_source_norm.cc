#include "calculator/cell/fission_source_norm.h"

#include "data/cross_sections.h"
#include "domain/finite_element_i.h"
#include "system/moments/spherical_harmonic_i.h"

namespace bart {

namespace calculator {

namespace cell {

template<int dim>
FissionSourceNorm<dim>::FissionSourceNorm(
    std::shared_ptr<domain::FiniteElementI<dim>> finite_element_ptr,
    std::shared_ptr<data::CrossSections> cross_sections_ptr)
    : finite_element_ptr_(finite_element_ptr),
      cross_sections_ptr_(cross_sections_ptr),
      cell_quadrature_points_(finite_element_ptr->n_cell_quad_pts()) {}


template<int dim>
double FissionSourceNorm<dim>::GetCellNorm(
    domain::CellPtr<dim> cell_ptr,
    system::moments::SphericalHarmonicI* system_moments) const {
  // TODO: change system_moments variable to system_moments_ptr
  const int material_id = cell_ptr->material_id();
  double fission_source = 0;

  if (cross_sections_ptr_->is_material_fissile.at(material_id)) {
    finite_element_ptr_->SetCell(cell_ptr);

    const int total_groups = system_moments->total_groups();
    const auto nu_sigma_f = cross_sections_ptr_->nu_sigma_f.at(material_id);

    for (int group = 0; group < total_groups; ++group) {

      auto scalar_flux_at_cell_quadrature =
          finite_element_ptr_->ValueAtQuadrature(
              system_moments->GetMoment({group, 0, 0}));

      for (int q = 0; q < cell_quadrature_points_; ++q) {
        double scalar_flux = scalar_flux_at_cell_quadrature.at(q) *
            finite_element_ptr_->Jacobian(q);
        fission_source += nu_sigma_f.at(group) * scalar_flux;
      }
    }
  }

  return fission_source;
}

template class FissionSourceNorm<1>;
template class FissionSourceNorm<2>;
template class FissionSourceNorm<3>;

} // namespace cell

} // namespace calculator

} // namespace bart