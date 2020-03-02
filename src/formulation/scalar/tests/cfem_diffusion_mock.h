#ifndef BART_SRC_FORMULATION_SCALAR_TESTS_CFEM_DIFFUSION_MOCK_H_
#define BART_SRC_FORMULATION_SCALAR_TESTS_CFEM_DIFFUSION_MOCK_H_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_accessor.h>

#include "system/moments/spherical_harmonic_types.h"
#include "formulation/scalar/diffusion_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace formulation {

namespace scalar {

template <int dim>
class CFEM_DiffusionMock : public DiffusionI<dim> {
 public:
  using typename DiffusionI<dim>::InitializationToken;
  using typename DiffusionI<dim>::BoundaryType;

  using typename DiffusionI<dim>::CellPtr;
  using typename DiffusionI<dim>::Matrix;
  using typename DiffusionI<dim>::Vector;
  using typename DiffusionI<dim>::GroupNumber;
  using typename DiffusionI<dim>::FaceNumber;

  MOCK_METHOD1_T(Precalculate, InitializationToken(const CellPtr& cell_ptr));
  MOCK_CONST_METHOD4_T(FillCellStreamingTerm, void(Matrix&,
      const InitializationToken,
      const CellPtr&,
      const GroupNumber));

  MOCK_CONST_METHOD4_T(FillCellCollisionTerm, void(Matrix&,
      const InitializationToken,
      const CellPtr&,
      const GroupNumber));

  MOCK_CONST_METHOD5_T(FillBoundaryTerm, void(Matrix&,
      const InitializationToken,
      const CellPtr&,
      const FaceNumber,
      const BoundaryType));

  MOCK_CONST_METHOD3_T(FillCellFixedSource, void(Vector& to_fill,
      const CellPtr&,
      const GroupNumber));

  MOCK_CONST_METHOD6_T(FillCellFissionSource, void(Vector&,
      const CellPtr&,
      const GroupNumber,
      const double,
      const system::moments::MomentVector&,
      const system::moments::MomentsMap&));

  MOCK_CONST_METHOD4_T(FillCellScatteringSource, void(Vector&,
      const CellPtr&,
      const GroupNumber,
      const system::moments::MomentsMap&));
};


} // namespace scalar

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_SCALAR_TESTS_CFEM_DIFFUSION_MOCK_H_