#include "calculator/cell/integrated_fission_source.hpp"

#include "data/cross_sections/material_cross_sections.hpp"
#include "domain/finite_element/finite_element_i.hpp"
#include "system/moments/spherical_harmonic_i.h"

namespace bart {

namespace calculator::cell {

namespace  {
template <int dim>
auto GetCellQuadraturePoints(domain::finite_element::FiniteElementI<dim>* finite_element_ptr) -> int {
  AssertThrow(finite_element_ptr != nullptr,
              dealii::ExcMessage("Error in constructor in IntegratedFissionSource: finite_element_ptr is null"))
  return finite_element_ptr->n_cell_quad_pts();
}
} // namespace

template<int dim>
IntegratedFissionSource<dim>::IntegratedFissionSource(std::shared_ptr<FiniteElement> finite_element_ptr,
                                                      std::shared_ptr<data::cross_sections::MaterialCrossSections> cross_sections_ptr)
    : finite_element_ptr_(finite_element_ptr),
      cross_sections_ptr_(cross_sections_ptr),
      cell_quadrature_points_(GetCellQuadraturePoints(finite_element_ptr.get())) {
  this->AssertPointerNotNull(cross_sections_ptr_.get(), "cross-sections", "IntegratedFissionSource constructor");
  this->set_description("integrated fission source calculator", utility::DefaultImplementation(true));
}

template<int dim>
auto IntegratedFissionSource<dim>::CellValue(CellPtr cell_ptr, MomentPtr system_moments_ptr) const -> double {
  const int material_id{ static_cast<int>(cell_ptr->material_id()) };
  double fission_source{ 0 };

  if (cross_sections_ptr_->is_material_fissile().at(material_id)) {
    finite_element_ptr_->SetCell(cell_ptr);

    const int total_groups = system_moments_ptr->total_groups();
    const auto nu_sigma_f = cross_sections_ptr_->nu_sigma_f().at(material_id);

    for (int group = 0; group < total_groups; ++group) {
      auto scalar_flux_at_cell_quadrature =
          finite_element_ptr_->ValueAtQuadrature(system_moments_ptr->GetMoment({group, 0, 0}));

      for (int q = 0; q < cell_quadrature_points_; ++q) {
        double scalar_flux = scalar_flux_at_cell_quadrature.at(q) * finite_element_ptr_->Jacobian(q);
        fission_source += nu_sigma_f.at(group) * scalar_flux;
      }
    }
  }

  return fission_source;
}

template class IntegratedFissionSource<1>;
template class IntegratedFissionSource<2>;
template class IntegratedFissionSource<3>;

} // namespace calculator::cell

} // namespace bart