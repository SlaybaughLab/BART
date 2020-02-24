#ifndef BART_SRC_FORMULATION_ANGULAR_TESTS_CFEM_SELF_ADJOINT_ANGULAR_FLUX_MOCK_H_
#define BART_SRC_FORMULATION_ANGULAR_TESTS_CFEM_SELF_ADJOINT_ANGULAR_FLUX_MOCK_H_

#include "formulation/angular/cfem_self_adjoint_angular_flux_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
class CFEMSelfAdjointAngularFluxMock :
    public CFEMSelfAdjointAngularFluxI<dim> {
 public:

  using typename CFEMSelfAdjointAngularFluxI<dim>::InitializationToken;
  MOCK_METHOD(void, FillBoundaryBilinearTerm, (FullMatrix& to_fill,
      const InitializationToken init_token,
      const domain::CellPtr<dim>& cell_ptr,
      const domain::FaceIndex,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number), (override));
  MOCK_METHOD(void, FillCellCollisionTerm, (FullMatrix&,
      const InitializationToken, const domain::CellPtr<dim>&,
      const system::EnergyGroup), (override));
  MOCK_METHOD(void, FillCellFissionSourceTerm, (Vector&,
      const InitializationToken, const domain::CellPtr<dim>&,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>>,
      const system::EnergyGroup, const double,
      const system::moments::MomentVector&,
      const system::moments::MomentsMap&), (override));
  MOCK_METHOD(void, FillCellFixedSourceTerm, (Vector&,
      const InitializationToken, const domain::CellPtr<dim>&,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>>,
      const system::EnergyGroup), (override));
  MOCK_METHOD(void, FillCellScatteringSourceTerm, (Vector&,
      const InitializationToken, const domain::CellPtr<dim>&,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>>,
      const system::EnergyGroup, const system::moments::MomentVector&,
      const system::moments::MomentsMap&), (override));
  MOCK_METHOD(void, FillCellStreamingTerm, (FullMatrix&,
      const InitializationToken, const domain::CellPtr<dim>&,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>>,
      const system::EnergyGroup), (override));
  MOCK_METHOD(InitializationToken, Initialize,
      (const domain::CellPtr<dim>&), (override));
};

} // namespace angular

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_TESTS_CFEM_SELF_ADJOINT_ANGULAR_FLUX_MOCK_H_
