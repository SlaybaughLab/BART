#ifndef BART_SRC_FORMULATION_SCALAR_TESTS_DRIFT_DIFFUSION_MOCK_HPP_
#define BART_SRC_FORMULATION_SCALAR_TESTS_DRIFT_DIFFUSION_MOCK_HPP_

#include "formulation/scalar/drift_diffusion_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::formulation::scalar {

template <int dim>
class DriftDiffusionMock : public DriftDiffusionI<dim> {
 public:
  using typename DriftDiffusionI<dim>::CellPtr;
  using typename DriftDiffusionI<dim>::EnergyGroup;
  using typename DriftDiffusionI<dim>::Matrix;
  using typename DriftDiffusionI<dim>::Vector;

  MOCK_METHOD(void, FillCellBoundaryTerm, (Matrix& to_fill, const CellPtr&,
      const domain::FaceIndex, const BoundaryType,
      std::function<Vector(const dealii::Tensor<1, dim>& normal_vector)> boundary_factor_function), (const, override));
  MOCK_METHOD(void, FillCellBoundaryTerm, (Matrix& to_fill, const CellPtr&, const domain::FaceIndex,
      const BoundaryType, const Vector& boundary_factor_at_global_dofs), (const, override));
  MOCK_METHOD(void, FillCellDriftDiffusionTerm, (Matrix& to_fill, const CellPtr&, const system::EnergyGroup,
                  const Vector& group_scalar_flux, (const std::array<Vector, dim>)& current), (const,override));
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_TESTS_DRIFT_DIFFUSION_MOCK_HPP_
