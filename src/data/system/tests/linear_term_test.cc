#include "data/system/term.h"

#include <memory>
#include <unordered_set>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/base/mpi.h>

#include "data/system/system_types.h"
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


} // namespace

