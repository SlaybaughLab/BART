#ifndef BART_SRC_FORMULATION_ANGULAR_SELF_ADJOINT_ANGULAR_FLUX_H_
#define BART_SRC_FORMULATION_ANGULAR_SELF_ADJOINT_ANGULAR_FLUX_H_

#include "data/cross_sections.h"
#include "domain/finite_element/finite_element_i.h"
#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "quadrature/quadrature_set_i.h"

#include <memory>

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
class SelfAdjointAngularFlux : public SelfAdjointAngularFluxI<dim> {
 public:

  SelfAdjointAngularFlux(
      std::shared_ptr<domain::finite_element::FiniteElementI<dim>>,
      std::shared_ptr<data::CrossSections>,
      std::shared_ptr<quadrature::QuadratureSetI<dim>>);

  void Initialize(const domain::CellPtr<dim>&) override;

  // Fill functions

  void FillBoundaryBilinearTerm(
      FullMatrix &to_fill,
      const domain::CellPtr<dim> &cell_ptr,
      const domain::FaceIndex face_number,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) override;

  void FillCellCollisionTerm(
      FullMatrix &to_fill,
      const domain::CellPtr<dim> &cell_ptr,
      const system::EnergyGroup group_number) override;

  void FillCellFissionSourceTerm(
      Vector &to_fill,
      const domain::CellPtr<dim> &cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const double k_eff,
      const system::moments::MomentVector &in_group_moment,
      const system::moments::MomentsMap &group_moments) override;

  void FillCellFixedSourceTerm(
      Vector &to_fill,
      const domain::CellPtr<dim> &cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) override;

  void FillCellScatteringSourceTerm(
      Vector &to_fill,
      const domain::CellPtr<dim> &cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const system::moments::MomentVector &in_group_moment,
      const system::moments::MomentsMap &group_moments) override;

  void FillCellStreamingTerm(
      FullMatrix &to_fill,
      const domain::CellPtr<dim> &cell_ptr,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) override;


  // Getters for pre-calculated values
  std::vector<double> OmegaDotGradient(int cell_quadrature_point,
                                       quadrature::QuadraturePointIndex) const;
  FullMatrix OmegaDotGradientSquared(int cell_quadrature_point,
                                     quadrature::QuadraturePointIndex) const;
  // Dependency getters
  domain::finite_element::FiniteElementI<dim>* finite_element_ptr() const {
    return finite_element_ptr_.get(); }
  data::CrossSections* cross_sections_ptr() const {
    return cross_sections_ptr_.get(); }
  quadrature::QuadratureSetI<dim>* quadrature_set_ptr() const {
    return quadrature_set_ptr_.get(); }

  std::map<int, FullMatrix> shape_squared() const {
    return shape_squared_; }

  bool is_initialized() const { return is_initialized_; }

 protected:
  // Validation Functions
  void ValidateMatrixSizeAndSetCell(const domain::CellPtr<dim>& cell_ptr,
                                    const FullMatrix& matrix_to_validate,
                                    std::string called_function_name) {
    ValidateAndSetCell(cell_ptr, called_function_name);
    ValidateMatrixSize(matrix_to_validate, called_function_name);
  }
  void ValidateVectorSizeAndSetCell(const domain::CellPtr<dim>& cell_ptr,
                                    const Vector& vector_to_validate,
                                    std::string called_function_name) {
    ValidateAndSetCell(cell_ptr, called_function_name);
    ValidateVectorSize(vector_to_validate, called_function_name);
  }
  void ValidateAndSetCell(const domain::CellPtr<dim>& cell_ptr,
                          std::string function_name);
  void ValidateMatrixSize(const FullMatrix&, std::string called_function_name);
  void ValidateVectorSize(const Vector&, std::string called_function_name);
  void VerifyInitialized(std::string called_function_name);

  // Combined implementation functions
  void FillCellSourceTerm(
      Vector& to_fill,
      const int material_id,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      std::vector<double> source);

  // Dependencies
  std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_ptr_;
  std::shared_ptr<data::CrossSections> cross_sections_ptr_;
  std::shared_ptr<quadrature::QuadratureSetI<dim>> quadrature_set_ptr_;
  // Geometric properties
  const int cell_degrees_of_freedom_ = 0; //!< Degrees of freedom per cell
  const int cell_quadrature_points_ = 0; //!< Quadrature points per cell
  const int face_quadrature_points_ = 0; //!< Quadrature points per face
  // Precalculated matrices and vectors
  using CellQuadratureIndex = int;
  using AngleIndex = int;
  std::map<std::pair<CellQuadratureIndex, AngleIndex>,
           dealii::Vector<double>> omega_dot_gradient_;
  std::map<std::pair<CellQuadratureIndex, AngleIndex>,
           FullMatrix> omega_dot_gradient_squared_;
  std::map<CellQuadratureIndex, FullMatrix> shape_squared_ = {};
  bool is_initialized_ = false;
};

} // namespace angular

} // namespace formulation

} //namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_SELF_ADJOINT_ANGULAR_FLUX_H_
