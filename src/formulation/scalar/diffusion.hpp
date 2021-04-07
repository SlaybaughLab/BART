#ifndef BART_SRC_FORMULATION_SCALAR_DIFFUSION_HPP_
#define BART_SRC_FORMULATION_SCALAR_DIFFUSION_HPP_

#include <memory>

#include <deal.II/lac/full_matrix.h>

#include "data/cross_sections/material_cross_sections.hpp"
#include "domain/finite_element/finite_element_i.hpp"
#include "formulation/scalar/diffusion_i.hpp"
#include "system/moments/spherical_harmonic_types.h"
#include "utility/has_dependencies.h"

namespace bart::formulation::scalar {

template <int dim>
class Diffusion : public DiffusionI<dim>, public utility::HasDependencies {
 public:
  using typename DiffusionI<dim>::BoundaryType;
  using typename DiffusionI<dim>::CellPtr;
  using typename DiffusionI<dim>::Matrix;
  using typename DiffusionI<dim>::Vector;
  using typename DiffusionI<dim>::GroupNumber;
  using typename DiffusionI<dim>::FaceNumber;

  // Dependency types
  using FiniteElement = typename domain::finite_element::FiniteElementI<dim>;
  using CrossSections = data::cross_sections::CrossSectionsI;

  // Moment types
  using MomentVector = system::moments::MomentVector;
  using MomentVectorMap = system::moments::MomentsMap;

  Diffusion(std::shared_ptr<FiniteElement>, std::shared_ptr<CrossSections>);

  auto Precalculate(const CellPtr& cell_ptr) -> void override;

  auto FillCellStreamingTerm(Matrix& to_fill, const CellPtr&, GroupNumber) const -> void override;

  auto FillCellCollisionTerm(Matrix& to_fill, const CellPtr&, GroupNumber) const -> void override;

  auto FillBoundaryTerm(Matrix& to_fill, const CellPtr&, FaceNumber, BoundaryType) const -> void override;

  auto FillCellFixedSource(Vector& to_fill, const CellPtr&, GroupNumber) const -> void override;

  auto FillCellFissionSource(Vector& to_fill, const CellPtr&, GroupNumber, double k_eigenvalue,
                             const MomentVector& in_group_moment,
                             const MomentVectorMap& group_moments) const -> void override;

  auto FillCellScatteringSource(Vector& to_fill, const CellPtr&, GroupNumber,
                                const MomentVectorMap & group_moments) const -> void override;

  auto is_initialized() const -> bool override { return is_initialized_; }

  // Getters & Setters
  /*! \brief Get precalculated matrices for the square of the shape function.
   *
   * \return Vector containing matrices corresponding to each quadrature point.
   */
  auto GetShapeSquared() const -> std::vector<Matrix> { return shape_squared_; }

  /*! \brief Get precalculated matrices for the square of the gradient
   * of the shape function.
   *
   * \return Vector containing matrices corresponding to each quadrature point.
   */
  auto GetGradientSquared() const -> std::vector<Matrix> { return gradient_squared_; }

 protected:
  //! Finite element object to provide shape function values
  std::shared_ptr<FiniteElement> finite_element_ptr_;
  //! Cross-sections object for cross-section data
  std::shared_ptr<CrossSections> cross_sections_ptr_;

  std::vector<Matrix> shape_squared_;
  std::vector<Matrix> gradient_squared_;
  const int cell_degrees_of_freedom_; //!< Number of degrees of freedom per cell
  const int cell_quadrature_points_; //!< Number of quadrature points per cell
  const int face_quadrature_points_; //!< Number of quadrature points per face
  bool is_initialized_{ false };

  auto VerifyInitialized(const std::string& called_function_name) const -> void;
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_DIFFUSION_HPP_