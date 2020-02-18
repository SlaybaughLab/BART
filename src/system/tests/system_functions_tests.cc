#include "system/system_functions.h"

#include "domain/tests/definition_mock.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "system/solution/mpi_group_angular_solution.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class SystemFunctionsTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(SystemFunctionsTests, bart::testing::AllDimensions);

TYPED_TEST(SystemFunctionsTests, SetUpMPIAngularSolutionTest) {
  constexpr int dim = this->dim;
  bart::system::solution::MPIGroupAngularSolutionMock mock_solution;
  domain::DefinitionMock<dim> mock_definition;

  bart::system::SetUpMPIAngularSolution<dim>(mock_solution,
                                             mock_definition);
}


} // namespace