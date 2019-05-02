#include "data/system/right_hand_side.h"

#include <memory>
#include <unordered_set>

#include <deal.II/base/mpi.h>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
using data::system::MPIVector;
using VariableTerms = data::system::RightHandSide::VariableTerms;

class SystemRightHandSideTest : public ::testing::Test {
 protected:
  SystemRightHandSideTest()
      : test_rhs({VariableTerms::kScatteringSource}) {};
  data::system::RightHandSide test_rhs;
  std::shared_ptr<MPIVector> test_vector, double_test_vector;
  void SetUp() override;
  void FillVector(MPIVector& to_fill, int multiple = 1);
};

void SystemRightHandSideTest::SetUp() {
  const int entries_per_processor = 10;
  const int n_processors =
      dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  test_vector = std::make_shared<MPIVector>(MPI_COMM_WORLD,
                                            entries_per_processor*n_processors,
                                            entries_per_processor);
  double_test_vector = std::make_shared<MPIVector>();
  double_test_vector->reinit(*test_vector);

  FillVector(*test_vector);
  FillVector(*double_test_vector, 2);
}

void SystemRightHandSideTest::FillVector(MPIVector& to_fill,
                                              const int multiple) {
  auto [first_entry, last_entry] = to_fill.local_range();
  for (unsigned int i = first_entry; i < last_entry; ++i) {
    to_fill(i) = i * multiple;
  }
  to_fill.compress(dealii::VectorOperation::insert);
}

TEST_F(SystemRightHandSideTest, Constructor) {
  std::unordered_set<VariableTerms> variable_terms =
      {VariableTerms::kScatteringSource};

  EXPECT_EQ(test_rhs.GetVariableTerms(), variable_terms);
}

TEST_F(SystemRightHandSideTest, SetFixedPtrTest) {
  test_rhs.SetFixedPtr({0,0}, test_vector);
  test_rhs.SetFixedPtr(1, double_test_vector);

  EXPECT_EQ(test_vector.use_count(), 2);
  EXPECT_EQ(double_test_vector.use_count(), 2);
}

TEST_F(SystemRightHandSideTest, GetFixedPtrIndexTest) {
  test_rhs.SetFixedPtr({0,0}, test_vector);
  test_rhs.SetFixedPtr({0,1}, double_test_vector);

  EXPECT_EQ(test_rhs.GetFixedPtr({0,0}), test_vector);
  EXPECT_EQ(test_rhs.GetFixedPtr({0,1}), double_test_vector);

  EXPECT_EQ(test_rhs.GetFixedPtr({2,0}), nullptr);
}

TEST_F(SystemRightHandSideTest, GetFixedPtrGroupTest) {
  test_rhs.SetFixedPtr(0, test_vector);
  test_rhs.SetFixedPtr({0,1}, double_test_vector);
  test_rhs.SetFixedPtr(1, double_test_vector);
  test_rhs.SetFixedPtr({1,1}, test_vector);

  EXPECT_EQ(test_rhs.GetFixedPtr(0), test_vector);
  EXPECT_EQ(test_rhs.GetFixedPtr({0,1}), double_test_vector);
  EXPECT_EQ(test_rhs.GetFixedPtr(1), double_test_vector);
  EXPECT_EQ(test_rhs.GetFixedPtr({1,1}), test_vector);

  EXPECT_EQ(test_rhs.GetFixedPtr(2), nullptr);
}

TEST_F(SystemRightHandSideTest, SetVariablePtrTest) {
  test_rhs.SetVariablePtr({0,0}, VariableTerms::kScatteringSource, test_vector);
  test_rhs.SetVariablePtr(1,
                          VariableTerms::kScatteringSource,
                          double_test_vector);
  EXPECT_EQ(test_vector.use_count(), 2);
  EXPECT_EQ(double_test_vector.use_count(), 2);

  EXPECT_ANY_THROW(test_rhs.SetVariablePtr({0,0},
                                           VariableTerms::kFissionSource,
                                           test_vector));
  EXPECT_ANY_THROW(test_rhs.SetVariablePtr(1,
                                           VariableTerms::kFissionSource,
                                           test_vector));
}

TEST_F(SystemRightHandSideTest, GetVariablePtrIndexTest) {
  auto term = VariableTerms::kScatteringSource;
  test_rhs.SetVariablePtr({0,0}, term, test_vector);
  test_rhs.SetVariablePtr({0,1}, term, double_test_vector);

  EXPECT_EQ(test_rhs.GetVariablePtr({0,0}, term), test_vector);
  EXPECT_EQ(test_rhs.GetVariablePtr({0,1}, term), double_test_vector);

  EXPECT_EQ(test_rhs.GetVariablePtr({2,0}, term), nullptr);
  EXPECT_ANY_THROW(test_rhs.GetVariablePtr({0,0},
                                           VariableTerms::kFissionSource));
}



} // namespace

