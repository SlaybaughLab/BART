#ifndef BART_SRC_FORMULATION_SCALAR_TWO_GRID_DIFFUSION_I_HPP_
#define BART_SRC_FORMULATION_SCALAR_TWO_GRID_DIFFUSION_I_HPP_

#include "data/cross_sections/one_group_cross_sections_i.hpp"
#include "formulation/scalar/diffusion.hpp"

namespace bart::formulation::scalar {

template <int dim>
class TwoGridDiffusion : public Diffusion<dim> {
 public:
  using typename Diffusion<dim>::CellPtr, typename Diffusion<dim>::Matrix, typename Diffusion<dim>::GroupNumber;
  using typename Diffusion<dim>::FaceNumber, typename Diffusion<dim>::BoundaryType;
  using CrossSections = data::cross_sections::CrossSectionsI;
  using FiniteElement = typename domain::finite_element::FiniteElementI<dim>;
  using MomentsMap = system::moments::MomentsMap;
  using MomentVector = system::moments::MomentVector;
  using OneGroupCrossSections = data::cross_sections::OneGroupCrossSectionsI;
  using Vector = dealii::Vector<double>;
  TwoGridDiffusion(std::shared_ptr<FiniteElement>, std::shared_ptr<CrossSections>,
                   std::shared_ptr<OneGroupCrossSections>);
  ~TwoGridDiffusion() = default;

  auto FillCellCollisionTerm(Matrix& to_fill, const CellPtr&, GroupNumber) const -> void override;
  auto FillCellFissionSource(Vector&, const CellPtr&, GroupNumber, double, const MomentVector&,
                             const MomentsMap&) const -> void override {};
  auto FillCellFixedSource(Vector&, const CellPtr&, GroupNumber) const -> void override {};
  auto FillCellScatteringSource(Vector&, const CellPtr&, GroupNumber, const MomentsMap&) const -> void override {};

  auto one_group_cross_sections_ptr() const { return one_group_cross_sections_ptr_.get(); }
 private:
  std::shared_ptr<OneGroupCrossSections> one_group_cross_sections_ptr_{ nullptr };
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_TWO_GRID_DIFFUSION_I_HPP_
