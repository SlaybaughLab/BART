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
  EXPECT_TRUE(test_convergence.isConverged(flux_one, flux_one));
}

TEST_F(FluxL1ThresholdTest, OneThresholdAway) {
  double to_add = flux_one.l1_norm() * 0.99 * test_convergence.GetThreshold();
  flux_two(2) += to_add;
  
  flux_two.compress(dealii::VectorOperation::values::add);
  EXPECT_TRUE(test_convergence.isConverged(flux_one, flux_two));
  EXPECT_TRUE(test_convergence.isConverged(flux_two, flux_one));
}

TEST_F(FluxL1ThresholdTest, TwoThresholdAway) {
  double to_add = flux_one.l1_norm() * 2*test_convergence.GetThreshold();
  flux_two(2) += to_add;
  flux_two.compress(dealii::VectorOperation::values::add);
  EXPECT_FALSE(test_convergence.isConverged(flux_one, flux_two));
  EXPECT_FALSE(test_convergence.isConverged(flux_two, flux_one));
}

TEST_F(FluxL1ThresholdTest, SetThreshold) {
  double to_set = 1e-5;
  test_convergence.SetThreshold(to_set);
  EXPECT_EQ(test_convergence.GetThreshold(), to_set);

  double to_add = flux_one.l1_norm() * 0.99 * to_set;
  flux_two(2) += to_add;
  
  flux_two.compress(dealii::VectorOperation::values::add);
  EXPECT_TRUE(test_convergence.isConverged(flux_one, flux_two));
  EXPECT_TRUE(test_convergence.isConverged(flux_two, flux_one));
  
}

