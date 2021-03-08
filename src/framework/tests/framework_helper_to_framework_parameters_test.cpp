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
  using K_EffectiveUpdaterName = eigenvalue::k_eigenvalue::K_EffectiveUpdaterName;

  // Test object
  std::unique_ptr<FrameworkHelper> test_helper_ptr_{ nullptr };

  // Supporting objects and mocks
  framework::FrameworkParameters default_parameters_;
  problem::ParametersMock parameters_mock_;
  SystemHelperMock* system_helper_obs_ptr_ { nullptr };

  // Test parameters
  const std::string material_mapping_filename_{ "test_data/mapping.test_material_mapping" };
  const std::string parsed_material_mapping_{ "1 1 2 2" };
  const std::unordered_map<int, std::string> material_filenames_{
    {0, "test_data/material/serialized/uo2_20.material"}};
  const int energy_groups_{ 7 }; // required by the material file
  const bool use_nda_{ false };

  auto SetExpectations(const framework::FrameworkParameters&, const std::string& material_mapping_filename) -> void;
  auto SetExpectations(const framework::FrameworkParameters&) -> void;
  auto SetUp() -> void override;
};

std::vector<int> ToIntVector(const std::vector<double>& to_convert) {
  return std::vector<int>(to_convert.begin(), to_convert.end());
}

auto FrameworkHelperToFrameworkParametersTest::SetUp() -> void {
  default_parameters_.output_filename_base = "test_output_name_" + std::to_string(test_helpers::RandomInt(1, 10));
  default_parameters_.neutron_energy_groups = energy_groups_;
  default_parameters_.equation_type = problem::EquationType::kDiffusion;
  default_parameters_.reflective_boundaries = {problem::Boundary::kYMax};
  default_parameters_.material_mapping = parsed_material_mapping_;
  default_parameters_.eigen_solver_type = problem::EigenSolverType::kPowerIteration;
  default_parameters_.k_effective_updater = K_EffectiveUpdaterName::kCalculatorViaFissionSource;
  default_parameters_.group_solver_type = problem::InGroupSolverType::kSourceIteration;
  default_parameters_.angular_quadrature_type = problem::AngularQuadType::kGaussLegendre;
  default_parameters_.angular_quadrature_order = quadrature::Order(test_helpers::RandomInt(1, 5));
  const auto dim = test_helpers::RandomInt(1, 3);
  using SpatialDimension = framework::FrameworkParameters::SpatialDimension;
  using DomainSize = framework::FrameworkParameters::DomainSize;
  using NumberOfCells = framework::FrameworkParameters::NumberOfCells;
  default_parameters_.spatial_dimension = SpatialDimension(dim);
  default_parameters_.domain_size = DomainSize(test_helpers::RandomVector(dim, 10.0, 20.0));
  default_parameters_.number_of_cells = NumberOfCells(ToIntVector(test_helpers::RandomVector(dim, 5, 10)));
  default_parameters_.uniform_refinements = test_helpers::RandomInt(1, 3);
  default_parameters_.discretization_type = problem::DiscretizationType::kContinuousFEM;
  default_parameters_.cell_finite_element_type = problem::CellFiniteElementType::kGaussian;
  using PolynomialDegree = framework::FrameworkParameters::PolynomialDegree;
  default_parameters_.polynomial_degree = PolynomialDegree(test_helpers::RandomInt(2, 5));
  default_parameters_.use_nda_ = use_nda_;

  test_helper_ptr_ = std::make_unique<FrameworkHelper>(std::make_unique<SystemHelperMock>());
  system_helper_obs_ptr_ = dynamic_cast<SystemHelperMock*>(test_helper_ptr_->system_helper_ptr());
}

auto FrameworkHelperToFrameworkParametersTest::SetExpectations(
    const framework::FrameworkParameters& parameters) -> void {
  SetExpectations(parameters, material_mapping_filename_);
}

auto FrameworkHelperToFrameworkParametersTest::SetExpectations(
    const framework::FrameworkParameters& parameters,
    const std::string& material_mapping_filename) -> void {
  using Boundary = problem::Boundary;
  EXPECT_CALL(parameters_mock_, OutputFilenameBase()).WillOnce(Return(parameters.output_filename_base));
  EXPECT_CALL(parameters_mock_, NEnergyGroups()).WillOnce(Return(parameters.neutron_energy_groups));
  EXPECT_CALL(parameters_mock_, TransportModel()).WillOnce(Return(parameters.equation_type));
  std::map<Boundary, bool> reflective_boundaries {
      {Boundary::kXMin, false}, {Boundary::kXMax, false}, {Boundary::kYMin, false}, {Boundary::kYMax, false},
      {Boundary::kZMin, false}, {Boundary::kZMax, false}};
  for (auto boundary : parameters.reflective_boundaries)
    reflective_boundaries.at(boundary) = true;
  EXPECT_CALL(parameters_mock_, ReflectiveBoundary()).WillOnce(Return(reflective_boundaries));
  EXPECT_CALL(parameters_mock_, MaterialMapFilename()).WillOnce(Return(material_mapping_filename));
  if (parameters.eigen_solver_type.has_value()) {
    EXPECT_CALL(parameters_mock_, EigenSolver()).WillOnce(Return(parameters.eigen_solver_type.value()));
  } else {
    EXPECT_CALL(parameters_mock_, EigenSolver()).WillOnce(Return(problem::EigenSolverType::kNone));
  }
  EXPECT_CALL(parameters_mock_, InGroupSolver()).WillOnce(Return(parameters.group_solver_type));
  EXPECT_CALL(parameters_mock_, AngularQuad()).WillOnce(Return(parameters.angular_quadrature_type));
  EXPECT_CALL(parameters_mock_, AngularQuadOrder()).WillOnce(Return(parameters.angular_quadrature_order.value().get()));
  EXPECT_CALL(parameters_mock_, SpatialDimension()).WillOnce(Return(parameters.spatial_dimension.get()));
  EXPECT_CALL(parameters_mock_, SpatialMax()).WillOnce(Return(parameters.domain_size.get()));
  EXPECT_CALL(parameters_mock_, NCells()).WillOnce(Return(parameters.number_of_cells.get()));
  EXPECT_CALL(parameters_mock_, UniformRefinements()).WillOnce(Return(parameters.uniform_refinements));
  EXPECT_CALL(parameters_mock_, Discretization()).WillOnce(Return(parameters.discretization_type));
  EXPECT_CALL(parameters_mock_, FEPolynomialDegree()).WillOnce(Return(parameters.polynomial_degree.get()));
  EXPECT_CALL(parameters_mock_, MaterialFilenames()).WillOnce(Return(material_filenames_));
  EXPECT_CALL(parameters_mock_, NumberOfMaterials()).WillOnce(Return(static_cast<int>(material_filenames_.size())));
  EXPECT_CALL(parameters_mock_, K_EffectiveUpdaterType()).WillOnce(Return(parameters.k_effective_updater));
  EXPECT_CALL(parameters_mock_, DoNDA()).WillOnce(Return(parameters.use_nda_));
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
  } else if (lhs.k_effective_updater != rhs.k_effective_updater) {
    return AssertionFailure() << "K-effective updaters do not match";
  } else if (lhs.use_nda_ != rhs.use_nda_) {
    return AssertionFailure() << "use NDA flag do not match";
  }
  return AssertionSuccess();
}

TEST_F(FrameworkHelperToFrameworkParametersTest, DefaultParameters) {
  auto test_parameters{ default_parameters_ };

  SetExpectations(test_parameters);

  auto returned_parameters = test_helper_ptr_->ToFrameworkParameters(parameters_mock_);
  EXPECT_TRUE(AreEqual(test_parameters, returned_parameters));
  ASSERT_TRUE(returned_parameters.cross_sections_.has_value());
  EXPECT_NE(returned_parameters.cross_sections_.value(), nullptr);
}

TEST_F(FrameworkHelperToFrameworkParametersTest, UseNDATrue) {
  auto test_parameters{ default_parameters_ };
  test_parameters.use_nda_ = true;
  SetExpectations(test_parameters);
  auto returned_parameters = test_helper_ptr_->ToFrameworkParameters(parameters_mock_);
  EXPECT_TRUE(AreEqual(returned_parameters, test_parameters));
}

TEST_F(FrameworkHelperToFrameworkParametersTest, SAAFWithLevelSymmetric) {
  auto test_parameters{ default_parameters_ };
  test_parameters.equation_type = problem::EquationType::kSelfAdjointAngularFlux;
  test_parameters.angular_quadrature_type = problem::AngularQuadType::kLevelSymmetricGaussian;

  SetExpectations(test_parameters);
  auto returned_parameters = test_helper_ptr_->ToFrameworkParameters(parameters_mock_);
  EXPECT_TRUE(AreEqual(test_parameters, returned_parameters));
  ASSERT_TRUE(returned_parameters.cross_sections_.has_value());
  EXPECT_NE(returned_parameters.cross_sections_.value(), nullptr);
}

TEST_F(FrameworkHelperToFrameworkParametersTest, KEffectiveUpdaterRayleighQuotient) {
  auto test_parameters{ default_parameters_ };
  test_parameters.k_effective_updater = eigenvalue::k_eigenvalue::K_EffectiveUpdaterName::kCalculatorViaRayleighQuotient;

  SetExpectations(test_parameters);
  auto returned_parameters = test_helper_ptr_->ToFrameworkParameters(parameters_mock_);
  EXPECT_TRUE(AreEqual(test_parameters, returned_parameters));
  ASSERT_TRUE(returned_parameters.cross_sections_.has_value());
  EXPECT_NE(returned_parameters.cross_sections_.value(), nullptr);
}


TEST_F(FrameworkHelperToFrameworkParametersTest, DiscontinuousFEm) {
  auto test_parameters{ default_parameters_ };
  test_parameters.discretization_type = problem::DiscretizationType::kDiscontinuousFEM;

  SetExpectations(test_parameters);
  auto returned_parameters = test_helper_ptr_->ToFrameworkParameters(parameters_mock_);
  EXPECT_TRUE(AreEqual(test_parameters, returned_parameters));
  ASSERT_TRUE(returned_parameters.cross_sections_.has_value());
  EXPECT_NE(returned_parameters.cross_sections_.value(), nullptr);
}

TEST_F(FrameworkHelperToFrameworkParametersTest, FixedSourceSolve) {
  auto test_parameters{ default_parameters_ };
  test_parameters.eigen_solver_type = std::nullopt;

  SetExpectations(test_parameters);
  auto returned_parameters = test_helper_ptr_->ToFrameworkParameters(parameters_mock_);
  EXPECT_TRUE(AreEqual(test_parameters, returned_parameters));
  ASSERT_TRUE(returned_parameters.cross_sections_.has_value());
  EXPECT_NE(returned_parameters.cross_sections_.value(), nullptr);
}

// Calls that throw

TEST_F(FrameworkHelperToFrameworkParametersTest, BadMaterialMapping) {
  using Boundary = problem::Boundary;
  auto test_parameters{ default_parameters_ };

  EXPECT_CALL(parameters_mock_, OutputFilenameBase()).WillOnce(Return(test_parameters.output_filename_base));
  EXPECT_CALL(parameters_mock_, NEnergyGroups()).WillOnce(Return(test_parameters.neutron_energy_groups));
  EXPECT_CALL(parameters_mock_, TransportModel()).WillOnce(Return(test_parameters.equation_type));
  EXPECT_CALL(parameters_mock_, InGroupSolver()).WillOnce(Return(test_parameters.group_solver_type));
  EXPECT_CALL(parameters_mock_, AngularQuad()).WillOnce(Return(test_parameters.angular_quadrature_type));
  EXPECT_CALL(parameters_mock_, AngularQuadOrder()).WillOnce(Return(test_parameters.angular_quadrature_order.value().get()));
  EXPECT_CALL(parameters_mock_, SpatialDimension()).WillOnce(Return(test_parameters.spatial_dimension.get()));
  EXPECT_CALL(parameters_mock_, SpatialMax()).WillOnce(Return(test_parameters.domain_size.get()));
  EXPECT_CALL(parameters_mock_, NCells()).WillOnce(Return(test_parameters.number_of_cells.get()));
  EXPECT_CALL(parameters_mock_, UniformRefinements()).WillOnce(Return(test_parameters.uniform_refinements));
  EXPECT_CALL(parameters_mock_, Discretization()).WillOnce(Return(test_parameters.discretization_type));
  EXPECT_CALL(parameters_mock_, FEPolynomialDegree()).WillOnce(Return(test_parameters.polynomial_degree.get()));
  std::map<Boundary, bool> reflective_boundaries {
      {Boundary::kXMin, false}, {Boundary::kXMax, false}, {Boundary::kYMin, false}, {Boundary::kYMax, false},
      {Boundary::kZMin, false}, {Boundary::kZMax, false}};
  for (auto boundary : test_parameters.reflective_boundaries)
    reflective_boundaries.at(boundary) = true;
  EXPECT_CALL(parameters_mock_, ReflectiveBoundary()).WillOnce(Return(reflective_boundaries));
  EXPECT_CALL(parameters_mock_, MaterialMapFilename()).WillOnce(Return("bad_material_file"));
  EXPECT_CALL(parameters_mock_, K_EffectiveUpdaterType()).WillOnce(Return(test_parameters.k_effective_updater));
  EXPECT_CALL(parameters_mock_, DoNDA()).WillOnce(Return(test_parameters.use_nda_));

  EXPECT_ANY_THROW({
    auto returned_parameters = test_helper_ptr_->ToFrameworkParameters(parameters_mock_);
                   });
}



} // namespace