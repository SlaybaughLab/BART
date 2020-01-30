#ifndef BART_SRC_FORMULATION_TESTS_ANGULAR_STAMPER_MOCK_H_
#define BART_SRC_FORMULATION_TESTS_ANGULAR_STAMPER_MOCK_H_

#include "formulation/angular_stamper_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace formulation {

template <int dim>
class AngularStamperMock : public AngularStamperI<dim> {
 public:
  MOCK_METHOD(void, StampBoundaryBilinearTerm, (
      system::MPISparseMatrix& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number), (override));

  MOCK_METHOD(void, StampCollisionTerm, (system::MPISparseMatrix& to_stamp,
      const system::EnergyGroup group_number), (override));

  MOCK_METHOD(void, StampFissionSourceTerm, (
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const double k_eff,
      const system::moments::MomentVector &in_group_moment,
      const system::moments::MomentsMap &group_moments), (override));

  MOCK_METHOD(void, StampFixedSourceTerm, (
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number), (override));

  MOCK_METHOD(void, StampScatteringSourceTerm, (
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const system::moments::MomentVector &in_group_moment,
      const system::moments::MomentsMap &group_moments), (override));

  MOCK_METHOD(void, StampStreamingTerm, (
      system::MPISparseMatrix& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number), (override));
};

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_TESTS_ANGULAR_STAMPER_MOCK_H_
