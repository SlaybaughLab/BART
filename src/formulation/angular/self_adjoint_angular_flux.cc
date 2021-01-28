#include "formulation/angular/self_adjoint_angular_flux.h"

#include <algorithm>
#include <sstream>
#include <cstdlib>

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
SelfAdjointAngularFlux<dim>::SelfAdjointAngularFlux(
    std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_ptr,
    std::shared_ptr<data::cross_sections::CrossSectionsI> cross_sections_ptr,
    std::shared_ptr<quadrature::QuadratureSetI<dim>> quadrature_set_ptr)
    : finite_element_ptr_(finite_element_ptr),
      cross_sections_ptr_(cross_sections_ptr),
      quadrature_set_ptr_(quadrature_set_ptr),
      cell_degrees_of_freedom_(finite_element_ptr->dofs_per_cell()),
      cell_quadrature_points_(finite_element_ptr->n_cell_quad_pts()),
      face_quadrature_points_(finite_element_ptr->n_face_quad_pts()) {}

template<int dim>
void SelfAdjointAngularFlux<dim>::Initialize(const domain::CellPtr<dim> &cell_ptr) {
  AssertThrow(cell_ptr.state() == dealii::IteratorState::valid,
              dealii::ExcMessage("Error in SelfAdjointAngularFlux Initialize, "
                                 "cell pointer is invalid."))

  finite_element_ptr_->SetCell(cell_ptr);
  shape_squared_ = {};
  omega_dot_gradient_ = {};

  /* Precalculated values are held in maps that are indexed by the cell
   * quadrature point where they are valid */
  for (int cell_quad_index = 0; cell_quad_index < cell_quadrature_points_;
       ++cell_quad_index) {
    formulation::FullMatrix shape_squared(cell_degrees_of_freedom_,
                                          cell_degrees_of_freedom_);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        shape_squared(i, j) =
            finite_element_ptr_->ShapeValue(i, cell_quad_index) *
            finite_element_ptr_->ShapeValue(j, cell_quad_index);
      }
    }
    shape_squared_.insert_or_assign(cell_quad_index, shape_squared);

    /* Each cell quadrature point has a further mapping based on the angular
     * quadrature point. */
    for (int angle_index : quadrature_set_ptr_->quadrature_point_indices()) {
      auto quadrature_point_ptr = quadrature_set_ptr_->GetQuadraturePoint(
          quadrature::QuadraturePointIndex(angle_index));

      dealii::Vector<double> omega_dot_gradient_vector(cell_degrees_of_freedom_);
      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
        omega_dot_gradient_vector[i] =
            quadrature_point_ptr->cartesian_position_tensor() *
            finite_element_ptr_->ShapeGradient(i, cell_quad_index);
      }

      omega_dot_gradient_.insert_or_assign({cell_quad_index, angle_index},
                                           omega_dot_gradient_vector);

      FullMatrix omega_dot_gradient_squared(cell_degrees_of_freedom_,
                                            cell_degrees_of_freedom_);
      omega_dot_gradient_squared.outer_product(omega_dot_gradient_vector,
                                               omega_dot_gradient_vector);
      omega_dot_gradient_squared_.insert_or_assign(
          {cell_quad_index, angle_index},
          omega_dot_gradient_squared);
    }
  }
  is_initialized_ = true;
}

template<int dim>
void SelfAdjointAngularFlux<dim>::FillBoundaryBilinearTerm(
    FullMatrix &to_fill,
    const domain::CellPtr<dim> &cell_ptr,
    const domain::FaceIndex face_number,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const system::EnergyGroup /*group_number*/) {
  VerifyInitialized(__FUNCTION__);
  ValidateMatrixSize(to_fill, __FUNCTION__);
  AssertThrow(cell_ptr.state() == dealii::IteratorState::valid,
              dealii::ExcMessage("Bad cell given to FilLBoundaryBilinearTerm"))
  finite_element_ptr_->SetFace(cell_ptr, face_number);

  auto normal_vector = finite_element_ptr_->FaceNormal();
  auto omega = quadrature_point->cartesian_position_tensor();

  const double normal_dot_omega = normal_vector * omega;

  if (normal_dot_omega > 0) {
    for (int f_q = 0; f_q < face_quadrature_points_; ++f_q) {
      const double jacobian = finite_element_ptr_->FaceJacobian(f_q);
      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
        for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
          to_fill(i,j) += normal_dot_omega
              * finite_element_ptr_->FaceShapeValue(i, f_q)
              * finite_element_ptr_->FaceShapeValue(j, f_q)
              * jacobian;
        }
      }
    }
  }
}

template<int dim>
auto SelfAdjointAngularFlux<dim>::FillReflectiveBoundaryLinearTerm(
    Vector& to_fill,
    const domain::CellPtr<dim>& cell_ptr,
    domain::FaceIndex face_number,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const dealii::Vector<double>& incoming_flux) -> double {
  VerifyInitialized(__FUNCTION__);
  ValidateVectorSize(to_fill, __FUNCTION__);
  AssertThrow(cell_ptr.state() == dealii::IteratorState::valid,
              dealii::ExcMessage("Bad cell given to FillReflectiveBoundaryLinearTerm"))
  finite_element_ptr_->SetFace(cell_ptr, face_number);

  auto normal_vector = finite_element_ptr_->FaceNormal();
  auto omega = quadrature_point->cartesian_position_tensor();
  double total_value_added{ 0 };

  const double normal_dot_omega = normal_vector * omega;

  if (normal_dot_omega < 0) {
    const auto incoming_angular_flux = finite_element_ptr_->ValueAtFaceQuadrature(
        incoming_flux);
    for (int f_q = 0; f_q < face_quadrature_points_; ++f_q) {
      const double jacobian = finite_element_ptr_->FaceJacobian(f_q);
      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
        const double value_to_add = normal_dot_omega
            * finite_element_ptr_->FaceShapeValue(i, f_q)
            * incoming_angular_flux.at(f_q)
            * jacobian;
        to_fill(i) -= value_to_add;
        total_value_added += std::abs(value_to_add);
      }
    }
  }
  return total_value_added;
}

template<int dim>
void SelfAdjointAngularFlux<dim>::FillCellCollisionTerm(
    FullMatrix &to_fill,
    const domain::CellPtr<dim> &cell_ptr,
    const system::EnergyGroup group_number) {
  VerifyInitialized(__FUNCTION__);
  ValidateMatrixSizeAndSetCell(cell_ptr, to_fill, __FUNCTION__);

  const int material_id = cell_ptr->material_id();
  const double sigma_t =
      cross_sections_ptr_->sigma_t().at(material_id).at(group_number.get());

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_ptr_->Jacobian(q);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) += sigma_t * finite_element_ptr_->ShapeValue(i, q) *
            finite_element_ptr_->ShapeValue(j, q) * jacobian;
      }
    }
  }
}

template<int dim>
auto SelfAdjointAngularFlux<dim>::FillCellFissionSourceTerm(
    Vector &to_fill,
    const domain::CellPtr<dim> & cell_ptr,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const system::EnergyGroup group_number,
    const double k_eff,
    const system::moments::MomentVector & in_group_moment,
    const system::moments::MomentsMap & group_moments) -> double {
  VerifyInitialized(__FUNCTION__);
  ValidateVectorSizeAndSetCell(cell_ptr, to_fill, __FUNCTION__);

  const int material_id = cell_ptr->material_id();
  const int group = group_number.get();
  double total_value_added{ 0 };

  if (cross_sections_ptr_->is_material_fissile().at(material_id)) {

    /* The scattering source is determined as the common values in both of the
     * scattering source terms in SAAF, specifically scalar flux times the
     * scattering cross-section per steradian */

    std::vector<double> fission_source(cell_quadrature_points_);

    // Get the contribution from each group
    for (const auto &moment_pair : group_moments) {
      auto &[index, moment] = moment_pair;
      const auto &[group_in, harmonic_l, harmonic_m] = index;

      if ((harmonic_l == 0) && (harmonic_m == 0)) {
        std::vector<double> scalar_flux(cell_quadrature_points_);

        if (group_in == group) {
          scalar_flux = finite_element_ptr_->ValueAtQuadrature(in_group_moment);
        } else {
          scalar_flux = finite_element_ptr_->ValueAtQuadrature(moment);
        }

        const auto fission_xfer_per_ster =
            cross_sections_ptr_->fiss_transfer_per_ster().at(material_id)(group_in,
                                                                        group);

        for (int q = 0; q < cell_quadrature_points_; ++q) {
          fission_source.at(q) +=
              fission_xfer_per_ster * scalar_flux.at(q) / k_eff;
        }
      }
    }

    total_value_added += std::abs(FillCellSourceTerm(to_fill, material_id, quadrature_point, group_number,
                                                     fission_source));
  }
  return total_value_added;
}

template<int dim>
void SelfAdjointAngularFlux<dim>::FillCellFixedSourceTerm(
    Vector &to_fill,
    const domain::CellPtr<dim> &cell_ptr,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const system::EnergyGroup group_number) {
  VerifyInitialized(__FUNCTION__);
  ValidateVectorSizeAndSetCell(cell_ptr, to_fill, __FUNCTION__);
  double q_per_ster = 0;
  const int material_id = cell_ptr->material_id();
  try {
    q_per_ster =
        cross_sections_ptr_->q_per_ster().at(material_id).at(group_number.get());
  } catch (std::out_of_range&) {
    return;
  }

  std::vector<double> fixed_source(cell_degrees_of_freedom_);
  std::fill(fixed_source.begin(), fixed_source.end(), q_per_ster);

  FillCellSourceTerm(to_fill, material_id, quadrature_point, group_number,
                     fixed_source);
}

template<int dim>
auto SelfAdjointAngularFlux<dim>::FillCellScatteringSourceTerm(
    Vector &to_fill,
    const domain::CellPtr<dim> &cell_ptr,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const system::EnergyGroup group_number,
    const system::moments::MomentVector &in_group_moment,
    const system::moments::MomentsMap &group_moments) -> double {
  VerifyInitialized(__FUNCTION__);
  ValidateVectorSizeAndSetCell(cell_ptr, to_fill, __FUNCTION__);

  const int material_id = cell_ptr->material_id();
  const int group = group_number.get();

  /* The scattering source is determined as the common values in both of the
   * scattering source terms in SAAF, specifically scalar flux times the
   * scattering cross-section per steradian */

  std::vector<double> scattering_source(cell_quadrature_points_);

  // Get the contribution from each group
  for (const auto& moment_pair : group_moments) {
    auto &[index, moment] = moment_pair;
    const auto &[group_in, harmonic_l, harmonic_m] = index;

    if ((harmonic_l == 0) && (harmonic_m == 0)) {
      std::vector<double> scalar_flux(cell_quadrature_points_);

      if (group_in == group) {
        scalar_flux = finite_element_ptr_->ValueAtQuadrature(in_group_moment);
      } else {
        scalar_flux = finite_element_ptr_->ValueAtQuadrature(moment);
      }

      const auto sigma_s_per_ster =
          cross_sections_ptr_->sigma_s_per_ster().at(material_id)(group, group_in);

      for (int q = 0; q < cell_quadrature_points_; ++q){
        scattering_source.at(q) += sigma_s_per_ster * scalar_flux.at(q);
      }
    }
  }

  return std::abs(FillCellSourceTerm(to_fill, material_id, quadrature_point, group_number, scattering_source));
}

template<int dim>
void SelfAdjointAngularFlux<dim>::FillCellStreamingTerm(
    FullMatrix &to_fill,
    const domain::CellPtr<dim> &cell_ptr,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const system::EnergyGroup group_number) {
  VerifyInitialized(__FUNCTION__);
  ValidateMatrixSizeAndSetCell(cell_ptr, to_fill, __FUNCTION__);

  const int material_id = cell_ptr->material_id();
  const double inverse_sigma_t =
      cross_sections_ptr_->inverse_sigma_t().at(material_id).at(group_number.get());
  const int angle_index = quadrature_set_ptr_->GetQuadraturePointIndex(
      quadrature_point);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_ptr_->Jacobian(q);
    const FullMatrix omega_dot_gradient_squared = OmegaDotGradientSquared(
        q, quadrature::QuadraturePointIndex(angle_index));
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) +=
            inverse_sigma_t * omega_dot_gradient_squared(i,j) * jacobian;
      }
    }
  }
}

template <int dim>
std::vector<double> SelfAdjointAngularFlux<dim>::OmegaDotGradient(
    int cell_quadrature_point,
    quadrature::QuadraturePointIndex angular_index) const {
  std::vector<double> return_vector(cell_degrees_of_freedom_);
  std::pair<int, int> index{cell_quadrature_point, angular_index.get()};
  for (int i = 0; i < cell_degrees_of_freedom_; ++i)
    return_vector.at(i) = omega_dot_gradient_.at(index)[i];
  return return_vector;
}

template <int dim>
FullMatrix SelfAdjointAngularFlux<dim>::OmegaDotGradientSquared(
    int cell_quadrature_point,
    quadrature::QuadraturePointIndex angular_index) const {
  return omega_dot_gradient_squared_.at(
      {cell_quadrature_point, angular_index.get()});
}

// PRIVATE FUNCTIONS ===========================================================
template <int dim>
void SelfAdjointAngularFlux<dim>::ValidateAndSetCell(
    const bart::domain::CellPtr<dim> &cell_ptr,
    std::string function_name) {
  std::string error{"Error in SelfAdjointAngularFlux function " +
      function_name + ": passed cell pointer is invalid"};
  AssertThrow(cell_ptr.state() == dealii::IteratorState::valid,
              dealii::ExcMessage(error))
  finite_element_ptr_->SetCell(cell_ptr);
}

template <int dim>
void SelfAdjointAngularFlux<dim>::ValidateMatrixSize(
    const bart::formulation::FullMatrix& to_validate,
    std::string called_function_name) {
  auto [rows, cols] = std::pair{to_validate.n_rows(), to_validate.n_cols()};

  std::ostringstream error_string;
  error_string << "Error in SelfAdjointAngularFlux function "
               << called_function_name
               <<": passed matrix size is invalid, expected size ("
               << cell_degrees_of_freedom_ << ", " << cell_degrees_of_freedom_
               << "), actual size: (" << rows << ", " << cols << ")";

  AssertThrow((static_cast<int>(rows) == cell_degrees_of_freedom_) &&
      (static_cast<int>(cols) == cell_degrees_of_freedom_),
      dealii::ExcMessage(error_string.str()))
}

template <int dim>
void SelfAdjointAngularFlux<dim>::ValidateVectorSize(
    const bart::formulation::Vector& to_validate,
    std::string called_function_name) {

  int rows = to_validate.size();

  std::ostringstream error_string;
  error_string << "Error in SelfAdjointAngularFlux function "
               << called_function_name
               <<": passed vector size is invalid, expected size ("
               << cell_degrees_of_freedom_ << ", 1), actual size: (" << rows
               << ", 1)";

  AssertThrow((static_cast<int>(rows) == cell_degrees_of_freedom_),
              dealii::ExcMessage(error_string.str()))
}



template <int dim>
auto SelfAdjointAngularFlux<dim>::FillCellSourceTerm(
    bart::formulation::Vector &to_fill,
    const int material_id,
    const std::shared_ptr<bart::quadrature::QuadraturePointI<dim>> quadrature_point,
    const bart::system::EnergyGroup group_number,
    std::vector<double> source) -> double {
  double total_value_added{ 0 };
  const double inverse_sigma_t =
      cross_sections_ptr_->inverse_sigma_t().at(material_id).at(group_number.get());
  const int angle_index = quadrature_set_ptr_->GetQuadraturePointIndex(
      quadrature_point);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_ptr_->Jacobian(q);
    const auto omega_dot_gradient = OmegaDotGradient(
        q, quadrature::QuadraturePointIndex(angle_index));

    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      const double value_to_add{ jacobian * source.at(q) * (finite_element_ptr_->ShapeValue(i, q) +
          omega_dot_gradient.at(i) * inverse_sigma_t)};
      to_fill(i) += value_to_add;
      total_value_added += std::abs(value_to_add);
    }
  }
  return total_value_added;
}

template<int dim>
void SelfAdjointAngularFlux<dim>::VerifyInitialized(
    std::string called_function_name) {
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

template class SelfAdjointAngularFlux<1>;
template class SelfAdjointAngularFlux<2>;
template class SelfAdjointAngularFlux<3>;

} // namespace angular

} // namespace formulation

} // namespace bart
