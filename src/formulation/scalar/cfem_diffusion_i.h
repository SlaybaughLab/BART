#ifndef BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_I_H_
#define BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_I_H_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_accessor.h>

#include "system/moments/spherical_harmonic_types.h"
#include "formulation/scalar/cfem_i.h"

namespace bart {

namespace formulation {

namespace scalar {

template <int dim>
class CFEM_DiffusionI : public CFEM_I {
 public:
  struct InitializationToken{};

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


  virtual ~CFEM_DiffusionI() = default;

  virtual InitializationToken Precalculate(const CellPtr& cell_ptr) = 0;

  virtual void FillCellStreamingTerm(Matrix& to_fill,
                             const InitializationToken,
                             const CellPtr& cell_ptr,
                             const GroupNumber group) const = 0;

  virtual void FillCellCollisionTerm(Matrix& to_fill,
                             const InitializationToken,
                             const CellPtr& cell_ptr,
                             const GroupNumber group) const = 0;

  virtual void FillBoundaryTerm(Matrix& to_fill,
                        const InitializationToken,
                        const CellPtr& cell_ptr,
                        const FaceNumber face_number,
                        const BoundaryType boundary_type) const = 0;

  virtual void FillCellFixedSource(Vector& to_fill,
                           const CellPtr& cell_ptr,
                           const GroupNumber group) const = 0;

  virtual void FillCellFissionSource(Vector& to_fill,
                             const CellPtr& cell_ptr,
                             const GroupNumber group,
                             const double k_effective,
                             const system::moments::MomentVector& in_group_moment,
                             const system::moments::MomentsMap& group_moments) const = 0;

  virtual void FillCellScatteringSource(Vector& to_fill,
                                const CellPtr& cell_ptr,
                                const GroupNumber group,
                                const system::moments::MomentsMap& group_moments) const = 0;

};

} // namespace scalar

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_SCALAR_CFEM_DIFFUSION_I_H_