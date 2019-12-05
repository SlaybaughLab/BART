#ifndef BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_
#define BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_

#include "data/cross_sections.h"
#include "domain/finite_element/finite_element_i.h"
#include "formulation/angular/cfem_self_adjoint_angular_flux_i.h"
#include "quadrature/quadrature_set_i.h"

#include <memory>

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
class CFEMSelfAdjointAngularFlux : public CFEMSelfAdjointAngularFluxI<dim> {
 public:
  using typename CFEMSelfAdjointAngularFluxI<dim>::InitializationToken;

  CFEMSelfAdjointAngularFlux(
      std::shared_ptr<domain::finite_element::FiniteElementI<dim>>,
      std::shared_ptr<data::CrossSections>,
      std::shared_ptr<quadrature::QuadratureSetI<dim>>);

  InitializationToken Initialize(const formulation::CellPtr<dim>&) override;

  // Fill functions
  void FillCellCollisionTerm(
      FullMatrix &to_fill,
      const InitializationToken init_token,
      const CellPtr<dim> &cell_ptr,
      const system::EnergyGroup group_number) override;

  void FillCellStreamingTerm(
      FullMatrix &to_fill,
      const InitializationToken init_token,
      const CellPtr<dim> &cell_ptr,
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

 protected:
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

};

} // namespace angular

} // namespace formulation

} //namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_
