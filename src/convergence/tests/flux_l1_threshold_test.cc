#include "../flux_l1_threshold.h"

#include <deal.II/lac/petsc_parallel_vector.h>

#include "../../test_helpers/test_helper_functions.h"
#include "../../test_helpers/gmock_wrapper.h"

class FluxL1ThresholdTest : public ::testing::Test {
 protected:
  bart::convergence::FluxL1Threshold test_convergence;
  bart::data::Flux flux_one;
  bart::data::Flux flux_two;
  void SetUp() override;
  
};

void FluxL1ThresholdTest::SetUp() {
  flux_one.reinit(MPI_COMM_WORLD, 5, 5);
  flux_two.reinit(MPI_COMM_WORLD, 5, 5);
  auto random_vector = btest::RandomVector(5, 0, 2);
  for (unsigned int i = 0; i < flux_one.size(); ++i) {
    flux_one(i) = random_vector[i];
    flux_two(i) = random_vector[i];
  }
  flux_one.compress(dealii::VectorOperation::values::insert);
  flux_two.compress(dealii::VectorOperation::values::insert);
}

TEST_F(FluxL1ThresholdTest, SameVector) {
  ASSERT_TRUE(test_convergence.isConverged(flux_one, flux_one));
}

TEST_F(FluxL1ThresholdTest, ExactlyThresholdAway) {
  flux_two(2) += test_convergence.GetThreshold();
  ASSERT_TRUE(test_convergence.isConverged(flux_one, flux_two));
  ASSERT_TRUE(test_convergence.isConverged(flux_two, flux_one));
}

TEST_F(FluxL1ThresholdTest, TwoThresholdAway) {
  flux_two(2) += 2*test_convergence.GetThreshold();
  ASSERT_FALSE(test_convergence.isConverged(flux_one, flux_two));
  ASSERT_FALSE(test_convergence.isConverged(flux_two, flux_one));
}



