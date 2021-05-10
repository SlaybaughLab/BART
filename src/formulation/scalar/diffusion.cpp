#include "formulation/scalar/diffusion.hpp"

namespace bart::formulation::scalar {

template<int dim>
Diffusion<dim>::Diffusion(std::shared_ptr<FiniteElement> finite_element_ptr,
                          std::shared_ptr<CrossSections> cross_sections_ptr)
    : finite_element_ptr_(finite_element_ptr),
      cross_sections_ptr_(cross_sections_ptr),
      cell_degrees_of_freedom_(finite_element_ptr_->dofs_per_cell()),
      cell_quadrature_points_(finite_element_ptr_->n_cell_quad_pts()),
      face_quadrature_points_(finite_element_ptr_->n_face_quad_pts()) {
  this->set_description("Diffusion Formulation", utility::DefaultImplementation(true));
  this->AssertPointerNotNull(finite_element_ptr_.get(), "finite element", "diffusion formulation constructor");
  this->AssertPointerNotNull(cross_sections_ptr_.get(), "cross sections", "diffusion formulation constructor");
}

template <int dim>
auto Diffusion<dim>::Precalculate(const CellPtr& cell_ptr) -> void {

  finite_element_ptr_->SetCell(cell_ptr);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    Matrix gradient_squared(cell_degrees_of_freedom_, cell_degrees_of_freedom_);
    Matrix shape_squared(cell_degrees_of_freedom_, cell_degrees_of_freedom_);

    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      const double shape_value_at_i{ finite_element_ptr_->ShapeValue(i, q) };
      const auto gradient_at_i{ finite_element_ptr_->ShapeGradient(i, q) };
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        shape_squared(i, j) = shape_value_at_i * finite_element_ptr_->ShapeValue(j, q);
        gradient_squared(i, j) =  gradient_at_i * finite_element_ptr_->ShapeGradient(j, q);
      }
    }
    shape_squared_.push_back(shape_squared);
    gradient_squared_.push_back(gradient_squared);
  }
  is_initialized_ = true;
}

template<int dim>
void Diffusion<dim>::FillCellConstantTerm(Vector &to_fill,
                                          const CellPtr &cell_ptr,
                                          const Vector &constant_vector) const {
  finite_element_ptr_->SetCell(cell_ptr);
  const auto constant_vector_at_quadrature{ this->finite_element_ptr_->ValueAtQuadrature(constant_vector) };

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double constant{ finite_element_ptr_->Jacobian(q) * constant_vector_at_quadrature.at(q) };
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      to_fill(i) +=  constant * finite_element_ptr_->ShapeValue(i, q);
    }
  }
}

template <int dim>
auto Diffusion<dim>::FillCellStreamingTerm(Matrix& to_fill, const CellPtr& cell_ptr,
                                           const GroupNumber group) const -> void {
  VerifyInitialized(__FUNCTION__);
  finite_element_ptr_->SetCell(cell_ptr);
  const int material_id = cell_ptr->material_id();

  const double diffusion_coef{ cross_sections_ptr_->diffusion_coef().at(material_id)[group] };

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double constant{ diffusion_coef * finite_element_ptr_->Jacobian(q) };
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) += constant * gradient_squared_[q](i, j);
      }
    }
  }
}

template <int dim>
auto Diffusion<dim>::FillCellCollisionTerm(Matrix& to_fill, const CellPtr& cell_ptr,
                                           const GroupNumber group) const -> void {
  VerifyInitialized(__FUNCTION__);
  finite_element_ptr_->SetCell(cell_ptr);
  const int material_id = cell_ptr->material_id();

  const double sigma_t{ cross_sections_ptr_->sigma_t().at(material_id)[group] };
  const double sigma_s{ cross_sections_ptr_->sigma_s().at(material_id)(group, group) };
  const double sigma_r{ sigma_t - sigma_s };

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double constant{ sigma_r * finite_element_ptr_->Jacobian(q) };
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) += constant * shape_squared_[q](i, j);
      }
    }
  }
}

template <int dim>
auto Diffusion<dim>::FillBoundaryTerm(Matrix& to_fill, const CellPtr& cell_ptr, const FaceNumber face_number,
                                      const BoundaryType boundary_type) const -> void {
  VerifyInitialized(__FUNCTION__);
  if (boundary_type == BoundaryType::kVacuum) {
    finite_element_ptr_->SetFace(cell_ptr, domain::FaceIndex(face_number));

    for (int q = 0; q < face_quadrature_points_; ++q) {
      const double constant{ 0.5 * finite_element_ptr_->FaceJacobian(q) };
      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
        const double face_shape_at_i{ finite_element_ptr_->FaceShapeValue(i, q) };
        for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
          to_fill(i, j) += constant * face_shape_at_i * finite_element_ptr_->FaceShapeValue(j, q);
        }
      }
    }
  }
}

template <int dim>
auto Diffusion<dim>::FillCellFixedSource(Vector& to_fill, const CellPtr& cell_ptr,
                                         const GroupNumber group) const -> void {

  finite_element_ptr_->SetCell(cell_ptr);
  int material_id = cell_ptr->material_id();
  double q{ 0 };
  try {
    q = cross_sections_ptr_->q().at(material_id).at(group);
  } catch (std::exception&) {
    return;
  }
  std::vector<double> cell_fixed_source(cell_quadrature_points_, q);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_ptr_->Jacobian(q);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      to_fill[i] += finite_element_ptr_->ShapeValue(i, q) * jacobian;
    }
  }
}

template <int dim>
auto Diffusion<dim>::FillCellFissionSource(Vector& to_fill, const CellPtr& cell_ptr, const GroupNumber group,
                                           const double k_eigenvalue, const MomentVector& in_group_moment,
                                           const MomentVectorMap & group_moments) const -> void {

  int material_id = cell_ptr->material_id();
  if (cross_sections_ptr_->is_material_fissile().at(material_id)) {
    finite_element_ptr_->SetCell(cell_ptr);

    std::vector<double> fission_source_at_quad_points(cell_quadrature_points_);

    // Get fission source contribution from each group at each quadrature point
    for (const auto& [index, moment] : group_moments) {
      const auto &[group_in, harmonic_l, harmonic_m] = index;
      if (harmonic_l == 0 && harmonic_m == 0) {
        std::vector<double> scalar_flux_at_quad_points(cell_quadrature_points_);

        if (group_in == group) {
          scalar_flux_at_quad_points = finite_element_ptr_->ValueAtQuadrature(in_group_moment);
        } else {
          scalar_flux_at_quad_points = finite_element_ptr_->ValueAtQuadrature(moment);
        }

        const auto fission_transfer{ cross_sections_ptr_->fiss_transfer().at(material_id)(group_in, group) };

        for (int q = 0; q < cell_quadrature_points_; ++q)
          fission_source_at_quad_points[q] += fission_transfer * scalar_flux_at_quad_points[q];
      }
    }

    // Integrate for each degree of freedom
    for (int q = 0; q < cell_quadrature_points_; ++q) {
      fission_source_at_quad_points[q] *= finite_element_ptr_->Jacobian(q) / k_eigenvalue;
      for (int i = 0; i < cell_degrees_of_freedom_; ++i)
        to_fill(i) += finite_element_ptr_->ShapeValue(i, q) * fission_source_at_quad_points[q];
    }
  }
}

template <int dim>
auto Diffusion<dim>::FillCellScatteringSource(Vector& to_fill, const CellPtr& cell_ptr, const GroupNumber group,
                                              const MomentVectorMap & group_moments) const -> void {
  finite_element_ptr_->SetCell(cell_ptr);
  const int material_id = cell_ptr->material_id();

  std::vector<double> scattering_source_at_quad_points(cell_quadrature_points_);

  // Get fission source contribution from each group at each quadrature point
  for (const auto& [index, moment] : group_moments) {
    const auto &[group_in, harmonic_l, harmonic_m] = index;

    // Check if scalar flux for an out-group
    if ((group_in != group) && (harmonic_l == 0) && (harmonic_m == 0)) {
      //std::vector<double> scalar_flux_at_quad_points(cell_quadrature_points_);

      const auto scalar_flux_at_quad_points{ finite_element_ptr_->ValueAtQuadrature(moment) };
      const double sigma_s{ cross_sections_ptr_->sigma_s().at(material_id)(group, group_in) };

      for (int q = 0; q < cell_quadrature_points_; ++q)
        scattering_source_at_quad_points[q] += sigma_s * scalar_flux_at_quad_points[q];
    }
  }

  // Integrate for each degree of freedom
  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double constant{ finite_element_ptr_->Jacobian(q) * scattering_source_at_quad_points.at(q) };
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      to_fill(i) += finite_element_ptr_->ShapeValue(i, q) * constant;
    }
  }
}

template<int dim>
auto Diffusion<dim>::VerifyInitialized(const std::string& called_function_name) const -> void {
  if (!is_initialized_) {
    std::ostringstream error_string;
    error_string << "Error in Diffusion function " << called_function_name << ": formulation has not been initialized, "
                 << "call Initialize prior to filing any cell matrices or vectors.";

    AssertThrow(false, dealii::ExcMessage(error_string.str()))
  }
}

template class Diffusion<1>;
template class Diffusion<2>;
template class Diffusion<3>;

} // namespace bart::formulation::scalar