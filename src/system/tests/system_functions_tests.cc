#include "system/system_functions.h"

#include "domain/tests/definition_mock.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "system/solution/mpi_group_angular_solution.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_assertions.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;
using ::testing::DoDefault;
using ::testing::Return, ::testing::ReturnRef;

void StampMPIVector(bart::system::MPIVector &to_fill, double value = 2) {
  auto [local_begin, local_end] = to_fill.local_range();
  for (unsigned int i = local_begin; i < local_end; ++i)
    to_fill(i) += value;
  to_fill.compress(dealii::VectorOperation::add);
}

template <typename DimensionWrapper>
class SystemFunctionsSetUpMPIAngularSolutionTests :
    public ::testing::Test,
    public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
  bart::system::solution::MPIGroupAngularSolutionMock mock_solution;
  domain::DefinitionMock<dim> mock_definition;

  const int total_angles_ = 3;
  std::map<bart::system::AngleIndex, bart::system::MPIVector> solution_map_;
  void SetUp() override;
};

template <typename DimensionWrapper>
void SystemFunctionsSetUpMPIAngularSolutionTests<DimensionWrapper>::SetUp() {
  this->SetUpDealii();
  ON_CALL(mock_solution, total_angles()).WillByDefault(Return(total_angles_));

  for (int i = 0; i < total_angles_; ++i) {
    bart::system::MPIVector mpi_vector;
    solution_map_.insert_or_assign(i, mpi_vector);
  }
  ON_CALL(mock_solution, solutions()).WillByDefault(ReturnRef(solution_map_));
  ON_CALL(mock_definition, locally_owned_dofs())
      .WillByDefault(Return(this->locally_owned_dofs_));
}

TYPED_TEST_SUITE(SystemFunctionsSetUpMPIAngularSolutionTests,
    bart::testing::AllDimensions);

TYPED_TEST(SystemFunctionsSetUpMPIAngularSolutionTests, BadNangles) {
  constexpr int dim = this->dim;

  std::array<int, 4> bad_total_angles{0, -1, 2, 4};

  for (const auto angle : bad_total_angles) {
    EXPECT_CALL(this->mock_solution, total_angles()).WillOnce(Return(angle));
    EXPECT_CALL(this->mock_solution, solutions()).WillOnce(DoDefault());
    EXPECT_ANY_THROW({
      bart::system::SetUpMPIAngularSolution<dim>(this->mock_solution,
                                                 this->mock_definition);
                     });
  }
}

TYPED_TEST(SystemFunctionsSetUpMPIAngularSolutionTests, SetUpDefaultValue) {
  constexpr int dim = this->dim;
  EXPECT_CALL(this->mock_solution, total_angles()).WillOnce(DoDefault());
  EXPECT_CALL(this->mock_solution, solutions()).WillOnce(DoDefault());
  EXPECT_CALL(this->mock_definition, locally_owned_dofs()).WillOnce(DoDefault());
  EXPECT_NO_THROW({
    bart::system::SetUpMPIAngularSolution<dim>(this->mock_solution,
                                               this->mock_definition);
  });
  bart::system::MPIVector expected_vector;
  expected_vector.reinit(this->locally_owned_dofs_, MPI_COMM_WORLD);
  StampMPIVector(expected_vector, 1.0);

  for (const auto& solution : this->solution_map_) {
    auto& mpi_vector = solution.second;
    ASSERT_GT(mpi_vector.size(), 0);
    EXPECT_TRUE(bart::testing::CompareMPIVectors(expected_vector, mpi_vector));
  }
}

TYPED_TEST(SystemFunctionsSetUpMPIAngularSolutionTests, ProvidedValues) {
  constexpr int dim = this->dim;
  const double value_to_set = btest::RandomDouble(0, 20);

  EXPECT_CALL(this->mock_solution, total_angles()).WillOnce(DoDefault());
  EXPECT_CALL(this->mock_solution, solutions()).WillOnce(DoDefault());
  EXPECT_CALL(this->mock_definition, locally_owned_dofs()).WillOnce(DoDefault());
  EXPECT_NO_THROW({
    bart::system::SetUpMPIAngularSolution<dim>(this->mock_solution,
                                               this->mock_definition,
                                               value_to_set);
                   });
  bart::system::MPIVector expected_vector;
  expected_vector.reinit(this->locally_owned_dofs_, MPI_COMM_WORLD);
  expected_vector = 0;
  StampMPIVector(expected_vector, value_to_set);

  for (const auto& solution : this->solution_map_) {
    auto& mpi_vector = solution.second;
    ASSERT_GT(mpi_vector.size(), 0);
    EXPECT_TRUE(bart::testing::CompareMPIVectors(expected_vector, mpi_vector));
  }
}




} // namespace