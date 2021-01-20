#include "formulation/scalar/drift_diffusion.hpp"
#include "formulation/scalar/scalar_formulation_factory.hpp"

namespace bart::formulation::scalar {

template<int dim>
DriftDiffusion<dim>::DriftDiffusion(std::shared_ptr<FiniteElement> finite_element_ptr,
                                    std::shared_ptr<CrossSections> cross_sections_ptr,
                                    std::shared_ptr<DriftDiffusionCalculator> drift_diffusion_calculator_ptr,
                                    std::shared_ptr<AngularFluxIntegrator> angular_flux_integrator_ptr)
    : finite_element_ptr_(finite_element_ptr),
      cross_sections_ptr_(cross_sections_ptr),
      drift_diffusion_calculator_ptr_(drift_diffusion_calculator_ptr),
      angular_flux_integrator_ptr_(angular_flux_integrator_ptr) {
  std::string function_name{"formulation::scalar::DriftDiffusion constructor"};
  AssertPointerNotNull(finite_element_ptr_.get(), "finite_element_ptr", function_name);
  AssertPointerNotNull(cross_sections_ptr_.get(), "cross_sections_ptr", function_name);
  AssertPointerNotNull(drift_diffusion_calculator_ptr_.get(), "drift_diffusion_calculator_ptr", function_name);
  AssertPointerNotNull(angular_flux_integrator_ptr_.get(), "angular flux integrator ptr", function_name);
  cell_quadrature_points_ = finite_element_ptr->n_cell_quad_pts();
  cell_degrees_of_freedom_ = finite_element_ptr->dofs_per_cell();
  face_quadrature_points_ = finite_element_ptr->n_face_quad_pts();
}

template <int dim>
bool DriftDiffusion<dim>::is_registered_ =
    DriftDiffusionIFactory<dim, std::shared_ptr<FiniteElement>, std::shared_ptr<CrossSections>,
                           std::shared_ptr<DriftDiffusionCalculator>, std::shared_ptr<AngularFluxIntegrator>>::get()
        .RegisterConstructor(DriftDiffusionFormulationName::kDefaultImplementation,
                             [](std::shared_ptr<FiniteElement> finite_element_ptr,
                                std::shared_ptr<CrossSections> cross_sections_ptr,
                                std::shared_ptr<DriftDiffusionCalculator> drift_diffusion_calculator_ptr,
                                std::shared_ptr<AngularFluxIntegrator> angular_flux_integrator_ptr)
                                -> std::unique_ptr<DriftDiffusionI<dim>> {
                               return std::make_unique<DriftDiffusion<dim>>(finite_element_ptr,
                                                                            cross_sections_ptr,
                                                                            drift_diffusion_calculator_ptr,
                                                                            angular_flux_integrator_ptr);
                             }
        );

template<int dim>
void DriftDiffusion<dim>::FillCellBoundaryTerm(Matrix& to_fill,
                                               const CellPtr& cell_ptr,
                                               domain::FaceIndex face_index,
                                               const BoundaryType boundary_type,
                                               const VectorMap& group_angular_flux) const {
  std::string error_prefix{"Error in DriftDiffusion<dim>::FillCellBoundaryTerm: "};
  AssertThrow(static_cast<int>(to_fill.m()) == cell_quadrature_points_, dealii::ExcMessage("matrix to fill has wrong m()"))
  AssertThrow(static_cast<int>(to_fill.n()) == cell_quadrature_points_, dealii::ExcMessage("matrix to fill has wrong n()"))

  if (boundary_type == BoundaryType::kVacuum) {
    finite_element_ptr_->SetFace(cell_ptr, face_index);
    auto normal_tensor = finite_element_ptr_->FaceNormal();
    dealii::Vector<double> normal_vector(dim);
    for (int dir = 0; dir < dim; ++dir)
      normal_vector[dir] = normal_tensor[dir];
    const auto directional_current = this->angular_flux_integrator_ptr_->DirectionalCurrent(group_angular_flux,
                                                                                            normal_vector);
    const auto directional_flux = this->angular_flux_integrator_ptr_->DirectionalFlux(group_angular_flux,
                                                                                      normal_vector);
    dealii::Vector<double> boundary_factor_at_global_dofs(directional_current.size());
    for (unsigned int i = 0; i < directional_current.size(); ++i) {
      const double flux{ directional_flux.at(i) };
      if (flux == 0) {
        boundary_factor_at_global_dofs[i] = 0;
      } else {
        boundary_factor_at_global_dofs[i] = directional_current.at(i)/directional_flux.at(i);
      }
    }

    auto boundary_factor_at_q{finite_element_ptr_->ValueAtFaceQuadrature(boundary_factor_at_global_dofs)};

    for (int face_q = 0; face_q < face_quadrature_points_; ++face_q) {
      const double jacobian{finite_element_ptr_->FaceJacobian(face_q)};
      for (int dof_i = 0; dof_i < cell_degrees_of_freedom_; ++dof_i) {
        auto shape_value_i{finite_element_ptr_->FaceShapeValue(dof_i, face_q)};
        for (int dof_j = 0; dof_j < cell_degrees_of_freedom_; ++dof_j) {
          to_fill(dof_i, dof_j) += shape_value_i * finite_element_ptr_->FaceShapeValue(dof_j, face_q)
              * boundary_factor_at_q.at(face_q) * jacobian;
        }
      }
    }
  }
}

template<int dim>
auto DriftDiffusion<dim>::FillCellDriftDiffusionTerm(Matrix &to_fill,
                                                     const CellPtr &cell_ptr,
                                                     system::EnergyGroup group,
                                                     const Vector &group_scalar_flux,
                                                     const std::array<Vector, dim> &current) const -> void{
  std::string error_prefix{"Error in DriftDiffusion<dim>::FillCellDriftDiffusionTerm: "};
  AssertThrow(static_cast<int>(to_fill.m()) == cell_quadrature_points_, dealii::ExcMessage("matrix to fill has wrong m()"))
  AssertThrow(static_cast<int>(to_fill.n()) == cell_quadrature_points_, dealii::ExcMessage("matrix to fill has wrong n()"))
  finite_element_ptr_->SetCell(cell_ptr);
  const auto material_id{ cell_ptr->material_id() };
  const double diffusion_coeff{ cross_sections_ptr_->diffusion_coef.at(material_id).at(group.get()) };
  const auto scalar_flux_at_q{ finite_element_ptr_->ValueAtQuadrature(group_scalar_flux) };

  std::array<std::vector<double>, dim> current_components_at_q;
  for (int i = 0; i < dim; ++i) {
    current_components_at_q.at(i) = finite_element_ptr_->ValueAtQuadrature(current.at(i));
  }

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian{ finite_element_ptr_->Jacobian(q) };
    dealii::Tensor<1, dim> current_vector_at_quadrature_point;
    for (int dir = 0; dir < dim; ++dir)
      current_vector_at_quadrature_point[dir] = current_components_at_q.at(dir).at(q);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      const auto drift_diffusion{drift_diffusion_calculator_ptr_->DriftDiffusionVector(
          scalar_flux_at_q.at(q),
          current_vector_at_quadrature_point,
          finite_element_ptr_->ShapeGradient(i, q),
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
