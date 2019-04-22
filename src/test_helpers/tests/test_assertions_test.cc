#include "test_helpers/test_assertions.h"

#include <vector>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using ::testing::AssertionResult;
using ::testing::AssertionFailure;
using ::testing::AssertionSuccess;

class TestAssertionsTest : public ::testing::Test {
 protected:
  std::vector<double> vector_1;
  std::vector<double> vector_2;
  dealii::Vector<double> dealii_vector_1;
  dealii::Vector<double> dealii_vector_2;

  void SetUp() override;
};

void TestAssertionsTest::SetUp() {
  vector_1 = btest::RandomVector(5, 0, 1.0);
  vector_2 = btest::RandomVector(5, 1.0, 2.0);
  dealii_vector_1.reinit(5);
  dealii_vector_2.reinit(5);

  for (int i = 0; i < 5; ++i) {
    dealii_vector_1[i] = vector_1[i];
    dealii_vector_2[i] = vector_2[i];
  }
}

TEST_F(TestAssertionsTest, GoodComparisonStdVector) {
  EXPECT_EQ(AssertionSuccess(),
            bart::testing::CompareVector(vector_1, vector_1));
}

TEST_F(TestAssertionsTest, BadComparisonStdVector) {
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareVector(vector_1, vector_2));
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareVector(vector_2, vector_1));
}

TEST_F(TestAssertionsTest, GoodComparisonDealiiVector) {
  EXPECT_EQ(AssertionSuccess(),
            bart::testing::CompareVector(dealii_vector_1, dealii_vector_1));
  EXPECT_EQ(AssertionSuccess(),
            bart::testing::CompareVector(dealii_vector_2, dealii_vector_2));
}

TEST_F(TestAssertionsTest, BadComparisonDealiiVector) {
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareVector(dealii_vector_1, dealii_vector_2));
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareVector(dealii_vector_2, dealii_vector_1));
}

} // namespace

