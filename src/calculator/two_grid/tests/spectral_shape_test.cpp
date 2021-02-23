#include "calculator/two_grid/spectral_shape.hpp"

#include <deal.II/lac/petsc_full_matrix.h>
#include <numeric>

#include "solver/eigenvalue/tests/spectral_radius_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;
using ::testing::Return, ::testing::ResultOf, ::testing::ContainerEq, ::testing::Pointwise, ::testing::DoubleNear;

class CalculatorTwoGridSpectralShapeTest : public ::testing::Test {
 public:
  using Matrix = dealii::FullMatrix<double>;
  using PETScFullMatrix = dealii::PETScWrappers::FullMatrix;
  using EigenvalueSolver = solver::eigenvalue::SpectralRadiusMock;
  using MatrixValues = std::vector<std::vector<double>>;
  using TestClass = calculator::two_grid::SpectralShape;

  // Test object
  std::unique_ptr<TestClass> test_class_{ nullptr };

  // Supporting objects
  Matrix sigma_t, sigma_s;
  PETScFullMatrix expected_a;

  // Supporting mocks and observation pointers
  EigenvalueSolver* eigenvalue_solver_obs_ptr_{ nullptr };


  // test parameters
  const MatrixValues expected_a_values{{0, 1.0/4.0, 1.0/2.0}, {0, 1.0/20.0, 7.0/10.0}, {0, 13.0/120.0, 31.0/60.0}};
  static constexpr int n_groups{ 3 };
  const double eigenvalue_{ test_helpers::RandomDouble(-100, 100) };
  std::vector<double> eigenvector_;
  std::vector<double> spectral_shape_function_;

  auto SetUp() -> void override;
};

auto CalculatorTwoGridSpectralShapeTest::SetUp() -> void {

  const MatrixValues sigma_t_values{{7, 0, 0}, {0, 10, 0}, {0, 0, 13}};
  const MatrixValues sigma_s_values{{3, 1, 2}, {1, 5, 3}, {2, 3, 7}};

  sigma_t.reinit(n_groups, n_groups);
  sigma_s.reinit(n_groups, n_groups);
  expected_a.reinit(n_groups, n_groups);

  for (int i = 0; i < n_groups; ++i) {
    for (int j = 0; j < n_groups; ++j) {
      sigma_t.set(i, j, sigma_t_values.at(i).at(j));
      sigma_s.set(i, j, sigma_s_values.at(i).at(j));
      expected_a.set(i, j, expected_a_values.at(i).at(j));
    }
  }
  // The eigenvector will be returned normalized in the L2 norm
  eigenvector_ = test_helpers::RandomVector(n_groups, 0, 100);

  double magnitude{ 0.0 };
  for (int i = 0; i < n_groups; ++i)
    magnitude += eigenvector_.at(i) * eigenvector_.at(i);

  std::for_each(eigenvector_.begin(), eigenvector_.end(), [magnitude](double& val) { val /= std::sqrt(magnitude); });
  double sum{ 0.0 };
  for (int i = 0; i < n_groups; ++i) {
    sum += eigenvector_.at(i);
  }
  // The spectral shape function should be normalized in the L1 norm
  spectral_shape_function_ = eigenvector_;
  for (int i = 0; i < n_groups; ++i) {
    spectral_shape_function_.at(i) /= sum;
  }

  auto eigenvalue_solver_ptr = std::make_unique<EigenvalueSolver>();
  eigenvalue_solver_obs_ptr_ = eigenvalue_solver_ptr.get();

  test_class_ = std::make_unique<TestClass>(std::move(eigenvalue_solver_ptr));
}

TEST_F(CalculatorTwoGridSpectralShapeTest, Getters) {
  ASSERT_NE(test_class_->eigenvalue_solver_ptr(), nullptr);
  EXPECT_EQ(test_class_->eigenvalue_solver_ptr(), eigenvalue_solver_obs_ptr_);
}

TEST_F(CalculatorTwoGridSpectralShapeTest, ConstructorThrowsOnNullDependency) {
  EXPECT_ANY_THROW(TestClass(nullptr););
}

auto vectorize_matrix(const dealii::PETScWrappers::MatrixBase& to_vectorize) {
  std::vector<std::vector<double>> return_vector;
  for (unsigned int i = 0; i < to_vectorize.m(); ++i) {
    std::vector<double> row_values;
    for (unsigned int j = 0; j < to_vectorize.n(); ++j) {
      row_values.push_back(to_vectorize(i, j));
    }
    return_vector.push_back(row_values);
  }
  return return_vector;
}

MATCHER(VectorsNear, "") {
  std::vector<double> result;
  auto v1 = std::get<0>(arg);
  auto v2 = std::get<1>(arg);
  std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(result), std::minus<double>());
  for (const auto val : result) {
    if (std::abs(val) > 1e-6)
      return false;
  }
  return true;
}

// SpectralShape should return the correct value
TEST_F(CalculatorTwoGridSpectralShapeTest, SpectralShape) {
  std::pair<double, std::vector<double>> eigenpair{eigenvalue_, eigenvector_};

  EXPECT_CALL(*this->eigenvalue_solver_obs_ptr_, SpectralRadius(ResultOf(vectorize_matrix,
                                                                         Pointwise(VectorsNear(), expected_a_values))))
      .WillOnce(Return(eigenpair));
  auto spectral_shape_vector = test_class_->CalculateSpectralShape(sigma_t, sigma_s);
  EXPECT_THAT(spectral_shape_vector, ContainerEq(spectral_shape_function_));
}

// Spectral shape should fix the returned vector if it is negative (it should be normalized already)
TEST_F(CalculatorTwoGridSpectralShapeTest, SpectralShapeNegative) {
  std::vector<double> negative_eigenvector(eigenvector_);
  std::for_each(negative_eigenvector.begin(), negative_eigenvector.end(), [](double& val) { val *= -1; });

  std::pair<double, std::vector<double>> eigenpair{eigenvalue_, negative_eigenvector};

  EXPECT_CALL(*this->eigenvalue_solver_obs_ptr_, SpectralRadius(ResultOf(vectorize_matrix,
                                                                         Pointwise(VectorsNear(), expected_a_values))))
      .WillOnce(Return(eigenpair));
  auto spectral_shape_vector = test_class_->CalculateSpectralShape(sigma_t, sigma_s);
  EXPECT_THAT(spectral_shape_vector, ContainerEq(spectral_shape_function_));
}

TEST_F(CalculatorTwoGridSpectralShapeTest, BadMatrixSize) {
  for (const int bad_size : {0, n_groups + 1, n_groups - 1} ) {
    Matrix bad_rows(bad_size, n_groups);
    Matrix bad_cols(n_groups, bad_size);
    for (auto matrix : {sigma_t, sigma_s}) {
      EXPECT_ANY_THROW(test_class_->CalculateSpectralShape(matrix, bad_rows));
      EXPECT_ANY_THROW(test_class_->CalculateSpectralShape(bad_rows, matrix));
      EXPECT_ANY_THROW(test_class_->CalculateSpectralShape(matrix, bad_cols));
      EXPECT_ANY_THROW(test_class_->CalculateSpectralShape(bad_cols, matrix));
    }
  }
}

} // namespace

