#include "formulation/scalar/diffusion.h"

namespace bart {

namespace formulation {

namespace scalar {

template<int dim>
Diffusion<dim>::Diffusion(std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element,
                          std::shared_ptr<data::cross_sections::CrossSections> cross_sections)
    : finite_element_(finite_element),
      cross_sections_(cross_sections),
      cell_degrees_of_freedom_(finite_element->dofs_per_cell()),
      cell_quadrature_points_(finite_element->n_cell_quad_pts()),
      face_quadrature_points_(finite_element->n_face_quad_pts()) {
  this->set_description("Diffusion Formulation",
                        utility::DefaultImplementation(true));
}

template <int dim>
void Diffusion<dim>::Precalculate(const CellPtr& cell_ptr) {

  finite_element_->SetCell(cell_ptr);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    Matrix gradient_squared(cell_degrees_of_freedom_,
                            cell_degrees_of_freedom_);
    Matrix shape_squared(cell_degrees_of_freedom_,
                         cell_degrees_of_freedom_);

    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        shape_squared(i, j) =
            finite_element_->ShapeValue(i, q) *
            finite_element_->ShapeValue(j, q);
        gradient_squared(i, j) =
            finite_element_->ShapeGradient(i, q) *
            finite_element_->ShapeGradient(j, q);
      }
    }
    shape_squared_.push_back(shape_squared);
    gradient_squared_.push_back(gradient_squared);
  }
  is_initialized_ = true;
}

template <int dim>
void Diffusion<dim>::FillCellStreamingTerm(Matrix& to_fill,
                                           const CellPtr& cell_ptr,
                                           const GroupNumber group) const {
  VerifyInitialized(__FUNCTION__);
  finite_element_->SetCell(cell_ptr);
  int material_id = cell_ptr->material_id();

  const double diffusion_coef =
      cross_sections_->diffusion_coef.at(material_id)[group];

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    double jacobian = finite_element_->Jacobian(q);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) += diffusion_coef * gradient_squared_[q](i, j) *jacobian;
      }
    }
  }

}

template <int dim>
void Diffusion<dim>::FillCellCollisionTerm(Matrix& to_fill,
                                           const CellPtr& cell_ptr,
                                           const GroupNumber group) const {
  VerifyInitialized(__FUNCTION__);
  finite_element_->SetCell(cell_ptr);
  int material_id = cell_ptr->material_id();

  const double sigma_t = cross_sections_->sigma_t.at(material_id)[group];
  const double sigma_s = cross_sections_->sigma_s.at(material_id)(group, group);
  double sigma_r = sigma_t - sigma_s;

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_->Jacobian(q);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) += sigma_r * shape_squared_[q](i, j) * jacobian;
      }
    }
  }
}

template <int dim>
void Diffusion<dim>::FillBoundaryTerm(Matrix& to_fill,
                                      const CellPtr& cell_ptr,
                                      const FaceNumber face_number,
                                      const BoundaryType boundary_type) const {
  VerifyInitialized(__FUNCTION__);
  if (boundary_type == BoundaryType::kVacuum) {
    finite_element_->SetFace(cell_ptr, domain::FaceIndex(face_number));

    for (int q = 0; q < face_quadrature_points_; ++q) {
      const double jacobian = finite_element_->FaceJacobian(q);
      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
        for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
          to_fill(i, j) +=
              0.5 * jacobian * finite_element_->FaceShapeValue(i, q)
                  * finite_element_->FaceShapeValue(j, q);
        }
      }
    }
  }
}

template <int dim>
void Diffusion<dim>::FillCellFixedSource(Vector& to_fill,
                                         const CellPtr& cell_ptr,
                                         const GroupNumber group) const {

  finite_element_->SetCell(cell_ptr);
  int material_id = cell_ptr->material_id();
  double q = 0;
  try {
    q = cross_sections_->q.at(material_id).at(group);
  } catch (std::exception&) {
    return;
  }
  std::vector<double> cell_fixed_source(cell_quadrature_points_, q);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_->Jacobian(q);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      to_fill[i] += finite_element_->ShapeValue(i, q) * jacobian;
    }
  }
}

template <int dim>
void Diffusion<dim>::FillCellFissionSource(
    Vector& to_fill,
    const CellPtr& cell_ptr,
    const GroupNumber group,
    const double k_effective,
    const system::moments::MomentVector& in_group_moment,
    const system::moments::MomentsMap& group_moments) const {

  int material_id = cell_ptr->material_id();
  if (cross_sections_->is_material_fissile.at(material_id)) {
    finite_element_->SetCell(cell_ptr);


    std::vector<double> fission_source_at_quad_points(cell_quadrature_points_);

    // Get fission source contribution from each group at each quadrature point
    for (const auto& moment_pair : group_moments) {
      auto &[index, moment] = moment_pair;
      int group_in = index[0];
      if (index[1] == 0 && index[2] == 0) {
        std::vector<double> scalar_flux_at_quad_points(cell_quadrature_points_);

        if (group_in == group) {
          scalar_flux_at_quad_points =
              finite_element_->ValueAtQuadrature(in_group_moment);
        } else {
          scalar_flux_at_quad_points =
              finite_element_->ValueAtQuadrature(moment);
        }

        auto fission_transfer =
            cross_sections_->fiss_transfer.at(material_id)(group_in, group);

        for (int q = 0; q < cell_quadrature_points_; ++q)
          fission_source_at_quad_points[q] +=
              fission_transfer * scalar_flux_at_quad_points[q];
      }
    }

    // Integrate for each degree of freedom
    for (int q = 0; q < cell_quadrature_points_; ++q) {
      fission_source_at_quad_points[q] *=
          finite_element_->Jacobian(q) / k_effective;

      for (int i = 0; i < cell_degrees_of_freedom_; ++i)
        to_fill(i) +=
            finite_element_->ShapeValue(i, q) * fission_source_at_quad_points[q];

    }
  }
}

template <int dim>
void Diffusion<dim>::FillCellScatteringSource(
      Vector& to_fill,
      const CellPtr& cell_ptr,
      const GroupNumber group,
      const system::moments::MomentsMap& group_moments) const {

  finite_element_->SetCell(cell_ptr);
  int material_id = cell_ptr->material_id();

  std::vector<double> scattering_source_at_quad_points(cell_quadrature_points_);

  // Get fission source contribution from each group at each quadrature point
  for (const auto& moment_pair : group_moments) {
    auto &[index, moment] = moment_pair;
    const auto &[group_in, harmonic_l, harmonic_m] = index;

    // Check if scalar flux for an out-group
    if ((group_in != group) && (harmonic_l == 0) && (harmonic_m == 0)) {
      std::vector<double> scalar_flux_at_quad_points(cell_quadrature_points_);

      scalar_flux_at_quad_points =
            finite_element_->ValueAtQuadrature(moment);

      const auto sigma_s =
          cross_sections_->sigma_s.at(material_id)(group, group_in);

      for (int q = 0; q < cell_quadrature_points_; ++q)
        scattering_source_at_quad_points[q] +=
            sigma_s * scalar_flux_at_quad_points[q];
    }
  }

  // Integrate for each degree of freedom
  for (int q = 0; q < cell_quadrature_points_; ++q) {
    scattering_source_at_quad_points[q] *= finite_element_->Jacobian(q);

    for (int i = 0; i < cell_degrees_of_freedom_; ++i)
      to_fill(i) +=
          finite_element_->ShapeValue(i, q) * scattering_source_at_quad_points[q];

  }
}

template<int dim>
void Diffusion<dim>::VerifyInitialized(std::string called_function_name) const {
  if (!is_initialized_) {

    std::ostringstream error_string;
    error_string << "Error in SelfAdjointAngularFlux function "
                 << called_function_name
                 << ": formulation has not been initialized, call Initialize "
                 << "prior to filing any cell matrices or vectors.";

    AssertThrow(false,
                dealii::ExcMessage(error_string.str()))
  }
}

template class Diffusion<1>;
template class Diffusion<2>;
template class Diffusion<3>;

} // namespace scalar

} // namespace formulation

} // namespace bart