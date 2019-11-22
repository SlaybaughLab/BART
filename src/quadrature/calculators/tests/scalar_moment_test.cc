#include "quadrature/calculators/scalar_moment.h"

#include <memory>

#include <deal.II/base/mpi.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include "system/system_types.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

using ::testing::Ref, ::testing::Return, ::testing::ReturnRef;

void SetVector(system::MPIVector& to_set, double value) {
  auto [first_row, last_row] = to_set.local_range();
  for (unsigned int i = first_row; i < last_row; ++i)
    to_set[i] = value;
  to_set.compress(dealii::VectorOperation::insert);
}

/* Tests for the SphericalHarmonicMomentsZerothMoment class. Mock quadrature set
 * is required.
 *
 * Test initial conditions: test object is constructed with mock quadrature set,
 * observation pointer is provided to set call expectations. Three mpi vectors
 * are provided with values 1, 10, and 100.
 */
class QuadratureCalculatorsScalarMomentTest : public ::testing::Test {
 protected:
  // Supporting objects
  system::solution::MPIGroupAngularSolutionMock mock_solution_;
  system::MPIVector mpi_vector_;


  // Test parameters
  const int n_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int n_entries_per_proc = 10;
  const double expected_result = 123.456;

  void SetUp() override;
};

void QuadratureCalculatorsScalarMomentTest::SetUp() {
  mpi_vector_.reinit(MPI_COMM_WORLD,
                     n_processes * n_entries_per_proc,
                     n_entries_per_proc);

  SetVector(mpi_vector_, 123.456);
}

TEST_F(QuadratureCalculatorsScalarMomentTest, BadNumberofAngles) {
  quadrature::calculators::ScalarMoment test_calculator;

  EXPECT_CALL(mock_solution_, total_angles())
      .WillOnce(Return(2));
  EXPECT_ANY_THROW(test_calculator.CalculateMoment(&mock_solution_, 0, 0, 0));
}

}