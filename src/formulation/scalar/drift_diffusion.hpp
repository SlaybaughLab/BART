#ifndef BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_HPP_
#define BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_HPP_

#include "formulation/scalar/drift_diffusion_i.hpp"

#include <memory>

#include "calculator/drift_diffusion/drift_diffusion_vector_calculator_i.hpp"
#include "data/cross_sections.h"
#include "domain/finite_element/finite_element_i.h"
#include "utility/has_dependencies.h"

namespace bart::formulation::scalar {

template <int dim>
class DriftDiffusion : public DriftDiffusionI<dim>, public utility::HasDependencies {
 public:
  using CrossSections = data::CrossSections;
  using DriftDiffusionCalculator = typename calculator::drift_diffusion::DriftDiffusionVectorCalculatorI<dim>;
  using FiniteElement = typename domain::finite_element::FiniteElementI<dim>;
  using typename DriftDiffusionI<dim>::Matrix;
  using typename DriftDiffusionI<dim>::CellPtr;
  using typename DriftDiffusionI<dim>::Vector;

  DriftDiffusion(std::shared_ptr<FiniteElement>,
                 std::shared_ptr<CrossSections>,
                 std::shared_ptr<DriftDiffusionCalculator>);
  [[deprecated]] auto FillCellDriftDiffusionTerm(Matrix& to_fill,
                                  const CellPtr& cell_ptr,
                                  system::EnergyGroup group,
                                  const Vector& group_scalar_flux,
                                  const Vector& integrated_angular_flux) const -> void override;

  auto FillCellDriftDiffusionTerm(Matrix &to_fill,
                                  const CellPtr &ptr,
                                  system::EnergyGroup group,
                                  const Vector &group_scalar_flux,
                                  const std::array<Vector, dim> &current) const -> void override;

  auto finite_element_ptr() const -> FiniteElement* { return finite_element_ptr_.get(); }
  auto cross_sections_ptr() const -> CrossSections* { return cross_sections_ptr_.get(); }
  auto drift_diffusion_calculator_ptr() const -> DriftDiffusionCalculator* {
    return drift_diffusion_calculator_ptr_.get(); }
 private:
  std::shared_ptr<FiniteElement> finite_element_ptr_{ nullptr };
  std::shared_ptr<CrossSections> cross_sections_ptr_{ nullptr };
  std::shared_ptr<DriftDiffusionCalculator> drift_diffusion_calculator_ptr_{ nullptr };
  int cell_quadrature_points_{ 0 };
  int cell_degrees_of_freedom_{ 0 };
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_HPP_
