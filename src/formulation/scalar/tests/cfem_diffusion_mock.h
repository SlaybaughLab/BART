#ifndef BART_SRC_FORMULATION_SCALAR_TESTS_CFEM_DIFFUSION_MOCK_H_
#define BART_SRC_FORMULATION_SCALAR_TESTS_CFEM_DIFFUSION_MOCK_H_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_accessor.h>

#include "data/moment_types.h"
#include "formulation/scalar/cfem_diffusion_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace formulation {

namespace scalar {

template <int dim>
class CFEM_DiffusionMock : public CFEM_DiffusionI<dim> {
 public:
  using typename CFEM_DiffusionI<dim>::InitializationToken;
  using typename CFEM_DiffusionI<dim>::BoundaryType;

  using typename CFEM_DiffusionI<dim>::CellPtr;
  using typename CFEM_DiffusionI<dim>::Matrix;
  using typename CFEM_DiffusionI<dim>::Vector;
  using typename CFEM_DiffusionI<dim>::GroupNumber;
  using typename CFEM_DiffusionI<dim>::FaceNumber;

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
      const data::MomentVector&,
      const data::MomentsMap&));

  MOCK_CONST_METHOD5_T(FillCellScatteringSource, void(Vector&,
      const CellPtr&,
      const GroupNumber,
      const data::MomentVector&,
      const data::MomentsMap&));
};


} // namespace scalar

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_SCALAR_TESTS_CFEM_DIFFUSION_MOCK_H_