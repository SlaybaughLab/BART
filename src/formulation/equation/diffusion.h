#ifndef BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
#define BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_

#include <memory>

#include "data/cross_sections.h"
#include "domain/finite_element_i.h"
#include "problem/parameter_types.h"
#include "formulation/types.h"
#include "formulation/equation/transport.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class Diffusion : public Transport<dim> {
 public:
  using MaterialID = int;
  using GroupNumber = int;
  using FaceNumber = int;
  using FullMatrix = dealii::FullMatrix<double>;
  using Vector = dealii::Vector<double>;
  using typename Transport<dim>::CellPtr;

  /*! \brief Constructor.
   *
   * Requires pointers to dependent objects required for assembly of the
   * equation.
   *
   * \param discretization type of discretization (CFEM/DFEM)
   * \param problem_type either a fixed source problem or an eigenvalue problem.
   * \param cross_sections object holding all material cross-sections.
   * \param finite_element object holding finite element information.
   */
  Diffusion(const DiscretizationType discretization,
            std::shared_ptr<data::CrossSections> cross_sections,
            std::shared_ptr<domain::FiniteElementI<dim>> finite_element)
      : Transport<dim>(EquationType::kScalar, discretization) {
    Transport<dim>::ProvideCrossSections(cross_sections);
    Transport<dim>::ProvideFiniteElement(finite_element);
  }

  virtual ~Diffusion() = default;

  void Precalculate(const CellPtr &cell_ptr);

  void FillCellStreamingTerm(FullMatrix &to_fill,
                             const CellPtr &cell_ptr,
                             const GroupNumber group) const;

  void FillCellCollisionTerm(FullMatrix &to_fill,
                             const CellPtr &cell_ptr,
                             const GroupNumber group) const;

  void FillBoundaryTerm(FullMatrix &to_fill,
                        const CellPtr &cell_ptr,
                        const GroupNumber group,
                        const FaceNumber face_number,
                        const BoundaryType boundary_type) const;

  void FillCellFixedSource(Vector &rhs_to_fill,
                           const CellPtr &cell_ptr,
                           const GroupNumber group) const;

  void FillCellScatteringSource(
      Vector &rhs_to_fill,
      const CellPtr &cell_ptr,
      const GroupNumber group,
      const data::ScalarFluxPtrs &scalar_flux_ptrs) const;

  void FillCellFissionSource(
      Vector &rhs_to_fill,
      const CellPtr &cell_ptr,
      const GroupNumber group,
      const double k_effective,
      const data::ScalarFluxPtrs &scalar_flux_ptrs) const;

 protected:
  std::vector<FullMatrix> shape_squared_;
  std::vector<FullMatrix> gradient_squared_;

  using Transport<dim>::SetCell;
  using Transport<dim>::SetFace;
  using Transport<dim>::finite_element_;
  using Transport<dim>::cross_sections_;
  using Transport<dim>::cell_degrees_of_freedom_;
  using Transport<dim>::cell_quadrature_points_;
  using Transport<dim>::face_quadrature_points_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_