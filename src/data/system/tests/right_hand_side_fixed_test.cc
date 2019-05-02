#include "data/system/right_hand_side_fixed.h"

#include <deal.II/base/mpi.h>
#include <memory>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
using data::system::MPIVector;

class SystemRightHandSideFixedTest : public ::testing::Test {
 protected:
  data::system::RightHandSideFixed test_rhs;
  std::shared_ptr<MPIVector> test_vector;
  void SetUp() override;
  void FillVector(MPIVector& to_fill, int multiple = 1);
};

void SystemRightHandSideFixedTest::SetUp() {
  const int entries_per_processor = 10;
  const int n_processors =
      dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  test_vector = std::make_shared<MPIVector>(MPI_COMM_WORLD,
                                            entries_per_processor*n_processors,
                                            entries_per_processor);
  FillVector(*test_vector);
}

void SystemRightHandSideFixedTest::FillVector(MPIVector& to_fill,
                                              const int multiple) {
  auto [first_entry, last_entry] = to_fill.local_range();
  for (unsigned int i = first_entry; i < last_entry; ++i) {
    to_fill(i) = i * multiple;
  }
  to_fill.compress(dealii::VectorOperation::insert);
}

TEST_F(SystemRightHandSideFixedTest, SetFixedPtrTest) {
  test_rhs.SetFixedPtr({0,0}, test_vector);

  EXPECT_EQ(test_vector.use_count(), 2);
}

TEST_F(SystemRightHandSideFixedTest, GetFixedPtrTest) {
  test_rhs.SetFixedPtr({0,0}, test_vector);

  EXPECT_EQ(test_vector, test_rhs.GetFixedPtr({0,0}));
  EXPECT_EQ(nullptr, test_rhs.GetFixedPtr({1,0}));
}

} // namespace

