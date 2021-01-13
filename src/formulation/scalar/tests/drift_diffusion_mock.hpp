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
  using typename DriftDiffusionI<dim>::VectorMap;

  MOCK_METHOD(void, FillCellBoundaryTerm, (Matrix& to_fill, const CellPtr&, const domain::FaceIndex,
      const BoundaryType, const VectorMap&), (const, override));
  MOCK_METHOD(void, FillCellDriftDiffusionTerm, (Matrix& to_fill, const CellPtr&, const system::EnergyGroup,
                  const Vector& group_scalar_flux, (const std::array<Vector, dim>)& current), (const,override));
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_TESTS_DRIFT_DIFFUSION_MOCK_HPP_
