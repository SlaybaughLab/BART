#include "framework/framework_helper.hpp"

#include "problem/tests/parameters_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "system/tests/system_helper_mock.hpp"

namespace  {

using namespace bart;
using ::testing::Return, ::testing::AssertionResult, ::testing::AssertionFailure, ::testing::AssertionSuccess;

class FrameworkHelperToFrameworkParametersTest : public ::testing::Test {
 public:
  static constexpr int dim{ 2 }; // dimension to run the tests in (should be arbitrary)
  using FrameworkHelper = typename framework::FrameworkHelper<dim>;
  using SystemHelperMock = const typename system::SystemHelperMock<dim>;

  // Test object
  std::unique_ptr<FrameworkHelper> test_helper_ptr_{ nullptr };

  // Supporting objects and mocks
  framework::FrameworkParameters default_parameters_;
  problem::ParametersMock parameters_mock_;
  SystemHelperMock* system_helper_obs_ptr_ { nullptr };

  auto SetExpectations(const framework::FrameworkParameters&) -> void;
  auto SetUp() -> void override;
};

auto FrameworkHelperToFrameworkParametersTest::SetUp() -> void {
  default_parameters_.output_filename_base = "test_output_name_" + std::to_string(test_helpers::RandomInt(1, 10));

  test_helper_ptr_ = std::make_unique<FrameworkHelper>(std::make_unique<SystemHelperMock>());
  system_helper_obs_ptr_ = dynamic_cast<SystemHelperMock*>(test_helper_ptr_->system_helper_ptr());
}

auto FrameworkHelperToFrameworkParametersTest::SetExpectations(
    const framework::FrameworkParameters& parameters) -> void {
  EXPECT_CALL(parameters_mock_, OutputFilenameBase()).WillOnce(Return(parameters.output_filename_base));
}

AssertionResult AreEqual(const framework::FrameworkParameters& lhs, const framework::FrameworkParameters& rhs) {
  if (lhs.name != rhs.name) {
    return AssertionFailure() << "names do not match";
  } else if (lhs.output_filename_base != rhs.output_filename_base) {
    return AssertionFailure() << "output filenames do not match";
  } else if (lhs.neutron_energy_groups != rhs.neutron_energy_groups) {
    return AssertionFailure() << "neutron energy groups do not match";
  } else if (lhs.equation_type != rhs.equation_type) {
    return AssertionFailure() << "equation types do not match";
  } else if (lhs.reflective_boundaries != rhs.reflective_boundaries) {
    return AssertionFailure() << "reflective boundaries do not match";
  } else if (lhs.material_mapping != rhs.material_mapping) {
    return AssertionFailure() << "material mappings do not match";
  } else if (lhs.eigen_solver_type != rhs.eigen_solver_type) {
    return AssertionFailure() << "eigen solver types do not match";
  } else if (lhs.group_solver_type != rhs.group_solver_type) {
    return AssertionFailure() << "group solver types do not match";
  } else if (lhs.angular_quadrature_type != rhs.angular_quadrature_type) {
    return AssertionFailure() << "angular quadrature types do not match";
  } else if (lhs.angular_quadrature_order != rhs.angular_quadrature_order) {
    return AssertionFailure() << "angular quadrature orders do not match";
  } else if (lhs.spatial_dimension != rhs.spatial_dimension) {
    return AssertionFailure() << "spatial dimension do not match";
  } else if (lhs.domain_size != rhs.domain_size) {
    return AssertionFailure() << "domain sizes do not match";
  } else if (lhs.number_of_cells != rhs.number_of_cells) {
    return AssertionFailure() << "number of cells do not match";
  } else if (lhs.uniform_refinements != rhs.uniform_refinements) {
    return AssertionFailure() << "uniform refinements do not match";
  } else if (lhs.discretization_type != rhs.discretization_type) {
    return AssertionFailure() << "discretization types do not match";
  } else if (lhs.cell_finite_element_type != rhs.cell_finite_element_type) {
    return AssertionFailure() << "cell finite element type do not match";
  } else if (lhs.polynomial_degree != rhs.polynomial_degree) {
    return AssertionFailure() << "polynomial degree do not match";
  }
  return AssertionSuccess();
}

TEST_F(FrameworkHelperToFrameworkParametersTest, Diffusion) {
  auto test_parameters{ default_parameters_ };

  SetExpectations(test_parameters);

  auto returned_parameters = test_helper_ptr_->ToFrameworkParameters(parameters_mock_);
  EXPECT_TRUE(AreEqual(test_parameters, returned_parameters));
}

} // namespace