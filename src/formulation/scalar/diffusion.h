#ifndef BART_SRC_FORMULATION_SCALAR_DIFFUSION_H_
#define BART_SRC_FORMULATION_SCALAR_DIFFUSION_H_

#include <memory>

#include <deal.II/lac/full_matrix.h>

#include "system/moments/spherical_harmonic_types.h"
#include "data/cross_sections.h"
#include "domain/finite_element/finite_element_i.h"
#include "formulation/scalar/diffusion_i.h"

namespace bart {

namespace formulation {

namespace scalar {

template <int dim>
class Diffusion : public DiffusionI<dim> {
 public:
  using typename DiffusionI<dim>::InitializationToken;
  using typename DiffusionI<dim>::BoundaryType;

  using typename DiffusionI<dim>::CellPtr;
  using typename DiffusionI<dim>::Matrix;
  using typename DiffusionI<dim>::Vector;
  using typename DiffusionI<dim>::GroupNumber;
  using typename DiffusionI<dim>::FaceNumber;


  Diffusion(std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element,
            std::shared_ptr<data::CrossSections> cross_sections);

  /*! \brief Precalculate matrices.
   *
   * \param cell_ptr any cell, no Jacobian is used so this is arbitrary.
   */
  InitializationToken Precalculate(const CellPtr& cell_ptr) override;

  void FillCellStreamingTerm(Matrix& to_fill,
                             const InitializationToken,
                             const CellPtr& cell_ptr,
                             const GroupNumber group) const override;

  void FillCellCollisionTerm(Matrix& to_fill,
                             const InitializationToken,
                             const CellPtr& cell_ptr,
                             const GroupNumber group) const override;

  void FillBoundaryTerm(Matrix& to_fill,
                        const InitializationToken,
                        const CellPtr& cell_ptr,
                        const FaceNumber face_number,
                        const BoundaryType boundary_type) const override;

  void FillCellFixedSource(Vector& to_fill,
                           const CellPtr& cell_ptr,
                           const GroupNumber group) const override;

  void FillCellFissionSource(Vector& to_fill,
                             const CellPtr& cell_ptr,
                             const GroupNumber group,
                             const double k_effective,
                             const system::moments::MomentVector& in_group_moment,
                             const system::moments::MomentsMap& group_moments) const override;

  void FillCellScatteringSource(Vector& to_fill,
                                const CellPtr& cell_ptr,
                                const GroupNumber group,
                                const system::moments::MomentsMap& group_moments) const override;

  // Getters & Setters
  /*! \brief Get precalculated matrices for the square of the shape function.
   *
   * \return Vector containing matrices corresponding to each quadrature point.
   */
  std::vector<Matrix> GetShapeSquared() const {
    return shape_squared_;
  }

  /*! \brief Get precalculated matrices for the square of the gradient
   * of the shape function.
   *
   * \return Vector containing matrices corresponding to each quadrature point.
   */
  std::vector<Matrix> GetGradientSquared() const {
    return gradient_squared_;
  }

 protected:
  //! Finite element object to provide shape function values
  std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_;
  //! Cross-sections object for cross-section data
  std::shared_ptr<data::CrossSections> cross_sections_;

  //Precalculated matrices
  std::vector<Matrix> shape_squared_;
  std::vector<Matrix> gradient_squared_;

  int cell_degrees_of_freedom_ = 0; //!< Number of degrees of freedom per cell
  int cell_quadrature_points_ = 0; //!< Number of quadrature points per cell
  int face_quadrature_points_ = 0; //!< Number of quadrature points per face
};


} // namespace scalar

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_SCALAR_DIFFUSION_H_