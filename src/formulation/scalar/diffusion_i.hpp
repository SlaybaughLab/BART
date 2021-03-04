#ifndef BART_SRC_FORMULATION_SCALAR_DIFFUSION_I_HPP_
#define BART_SRC_FORMULATION_SCALAR_DIFFUSION_I_HPP_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_accessor.h>

#include "system/moments/spherical_harmonic_types.h"
#include "utility/has_description.h"

//! Scalar (non-angular) formulations of the transport equation.
namespace bart::formulation::scalar {

template <int dim>
class DiffusionI : public utility::HasDescription {
 public:
  enum class BoundaryType {
    kVacuum,
    kReflective
  };

  //! Pointer to a cell iterator returned by a dof object.
  using CellPtr = typename dealii::DoFHandler<dim>::active_cell_iterator;
  using Matrix = dealii::FullMatrix<double>;
  using Vector = dealii::Vector<double>;
  using GroupNumber = int;
  using FaceNumber = int;

  virtual ~DiffusionI() = default;

  virtual auto Precalculate(const CellPtr& cell_ptr) -> void = 0;

  virtual auto FillCellStreamingTerm(Matrix& to_fill, const CellPtr&, GroupNumber) const -> void = 0;

  virtual auto FillCellCollisionTerm(Matrix& to_fill, const CellPtr&, GroupNumber) const -> void = 0;

  virtual auto FillBoundaryTerm(Matrix& to_fill, const CellPtr&, FaceNumber, BoundaryType) const -> void = 0;

  virtual auto FillCellFixedSource(Vector& to_fill, const CellPtr&, GroupNumber) const -> void = 0;

  virtual auto FillCellFissionSource(Vector& to_fill, const CellPtr&, GroupNumber, double k_eigenvalue,
                                     const system::moments::MomentVector& in_group_moment,
                                     const system::moments::MomentsMap& group_moments) const -> void = 0;

  virtual auto FillCellScatteringSource(Vector& to_fill, const CellPtr&, GroupNumber,
                                        const system::moments::MomentsMap& group_moments) const -> void = 0;

  virtual auto is_initialized() const -> bool = 0;
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_DIFFUSION_I_HPP_