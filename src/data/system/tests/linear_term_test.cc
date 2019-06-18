#include "data/system/term.h"

#include <memory>
#include <unordered_set>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/base/mpi.h>

#include "system/system_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
using data::system::MPIVector;
using VariableLinearTerms = data::system::VariableLinearTerms;

class SystemLinearTermTest : public ::testing::Test {
 protected:
  SystemLinearTermTest()
      : test_linear_term({VariableLinearTerms::kScatteringSource}) {};

  data::system::MPILinearTerm test_linear_term;

  std::shared_ptr<MPIVector> test_vector, double_test_vector;
  void SetUp() override;
  void FillVector(MPIVector& to_fill, int multiple = 1);
};

void SystemLinearTermTest::SetUp() {
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

void SystemLinearTermTest::FillVector(MPIVector& to_fill,
                                         const int multiple) {
  auto [first_entry, last_entry] = to_fill.local_range();
  for (unsigned int i = first_entry; i < last_entry; ++i) {
    to_fill(i) = i * multiple;
  }
  to_fill.compress(dealii::VectorOperation::insert);
}

TEST_F(SystemLinearTermTest, Constructor) {
  std::unordered_set<VariableLinearTerms> variable_terms =
      {VariableLinearTerms::kScatteringSource};

  EXPECT_EQ(test_linear_term.GetVariableTerms(), variable_terms);
}

TEST_F(SystemLinearTermTest, SetFixedPtrTest) {
  test_linear_term.SetFixedTermPtr({0, 0}, test_vector);
  test_linear_term.SetFixedTermPtr(1, double_test_vector);

  EXPECT_EQ(test_vector.use_count(), 2);
  EXPECT_EQ(double_test_vector.use_count(), 2);
}


TEST_F(SystemLinearTermTest, GetFixedPtrIndexTest) {
  test_linear_term.SetFixedTermPtr({0, 0}, test_vector);
  test_linear_term.SetFixedTermPtr({0, 1}, double_test_vector);

  EXPECT_EQ(test_linear_term.GetFixedTermPtr({0, 0}), test_vector);
  EXPECT_EQ(test_linear_term.GetFixedTermPtr({0, 1}), double_test_vector);

  EXPECT_EQ(test_linear_term.GetFixedTermPtr({2, 0}), nullptr);
}

TEST_F(SystemLinearTermTest, GetFixedPtrGroupTest) {
  test_linear_term.SetFixedTermPtr(0, test_vector);
  test_linear_term.SetFixedTermPtr({0, 1}, double_test_vector);
  test_linear_term.SetFixedTermPtr(1, double_test_vector);
  test_linear_term.SetFixedTermPtr({1, 1}, test_vector);

  EXPECT_EQ(test_linear_term.GetFixedTermPtr(0), test_vector);
  EXPECT_EQ(test_linear_term.GetFixedTermPtr({0, 1}), double_test_vector);
  EXPECT_EQ(test_linear_term.GetFixedTermPtr(1), double_test_vector);
  EXPECT_EQ(test_linear_term.GetFixedTermPtr({1, 1}), test_vector);

  EXPECT_EQ(test_linear_term.GetFixedTermPtr(2), nullptr);
}

TEST_F(SystemLinearTermTest, SetVariablePtrTest) {
  test_linear_term.SetVariableTermPtr({0, 0}, VariableLinearTerms::kScatteringSource, test_vector);
  test_linear_term.SetVariableTermPtr(1,
                              VariableLinearTerms::kScatteringSource,
                              double_test_vector);
  EXPECT_EQ(test_vector.use_count(), 2);
  EXPECT_EQ(double_test_vector.use_count(), 2);

  EXPECT_ANY_THROW(test_linear_term.SetVariableTermPtr({0, 0},
                                               VariableLinearTerms::kFissionSource,
                                               test_vector));
  EXPECT_ANY_THROW(test_linear_term.SetVariableTermPtr(1,
                                               VariableLinearTerms::kFissionSource,
                                               test_vector));
}

TEST_F(SystemLinearTermTest, GetVariablePtrIndexTest) {
  auto term = VariableLinearTerms::kScatteringSource;
  test_linear_term.SetVariableTermPtr({0, 0}, term, test_vector);
  test_linear_term.SetVariableTermPtr({0, 1}, term, double_test_vector);

  EXPECT_EQ(test_linear_term.GetVariableTermPtr({0, 0}, term), test_vector);
  EXPECT_EQ(test_linear_term.GetVariableTermPtr({0, 1}, term), double_test_vector);

  EXPECT_EQ(test_linear_term.GetVariableTermPtr({2, 0}, term), nullptr);
  EXPECT_ANY_THROW(test_linear_term.GetVariableTermPtr({0, 0},
                                               VariableLinearTerms::kFissionSource));
}

TEST_F(SystemLinearTermTest, GetVariablePtrGroupTest) {
  auto term = VariableLinearTerms::kScatteringSource;

  test_linear_term.SetVariableTermPtr(0, term, test_vector);
  test_linear_term.SetVariableTermPtr({0, 1}, term, double_test_vector);
  test_linear_term.SetVariableTermPtr(1, term, double_test_vector);
  test_linear_term.SetVariableTermPtr({1, 1}, term, test_vector);

  EXPECT_EQ(test_linear_term.GetVariableTermPtr(0, term), test_vector);
  EXPECT_EQ(test_linear_term.GetVariableTermPtr({0, 1}, term), double_test_vector);
  EXPECT_EQ(test_linear_term.GetVariableTermPtr(1, term), double_test_vector);
  EXPECT_EQ(test_linear_term.GetVariableTermPtr({1, 1}, term), test_vector);

  EXPECT_EQ(test_linear_term.GetVariableTermPtr(2, term), nullptr);
  EXPECT_ANY_THROW(test_linear_term.GetVariableTermPtr(0,
                                               VariableLinearTerms::kFissionSource));
}

} // namespace

