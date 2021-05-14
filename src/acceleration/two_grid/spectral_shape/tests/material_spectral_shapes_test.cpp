#include "acceleration/two_grid/spectral_shape/material_spectral_shapes.hpp"

#include <deal.II/lac/full_matrix.h>

#include "spectral_shape_mock.hpp"
#include "data/cross_sections/tests/cross_sections_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "material_spectral_shapes_mock.hpp"

namespace  {

using namespace bart;
using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::_;

class CalculatorTwoGridMaterialSpectralShapes : public ::testing::Test {
 public:
  using CrossSections = NiceMock<data::cross_sections::CrossSectionsMock>;
  using DealiiMatrix = dealii::FullMatrix<double>;
  using SpectralShapeMock = NiceMock<acceleration::two_grid::spectral_shape::SpectralShapeMock>;
  using MaterialSpectralShapes = acceleration::two_grid::spectral_shape::MaterialSpectralShapes;
  // Test object
  std::unique_ptr<MaterialSpectralShapes> test_calculator_{ nullptr };

  // Supporing mock objects and observation pointers
  SpectralShapeMock* spectral_shape_calculator_obs_ptr_{ nullptr };
  std::shared_ptr<CrossSections> cross_sections_ptr_{ std::make_shared<CrossSections>() };

  // Test parameters
  const int n_groups{ test_helpers::RandomInt(5, 10) };
  std::unordered_map<int, std::vector<double>> material_id_mapped_to_spectral_shape_;

  auto SetUp() -> void override;
  auto MakeMatrix(std::vector<double> values) const -> DealiiMatrix;
};

auto CalculatorTwoGridMaterialSpectralShapes::SetUp() -> void {
  const double max_cross_section_value{ 10.0 };
  const double max_spectral_shape_value{ 100 };
  auto spectral_shape_calculator_ptr = std::make_unique<SpectralShapeMock>();
  spectral_shape_calculator_obs_ptr_ = spectral_shape_calculator_ptr.get();

  std::unordered_map<int, DealiiMatrix> sigma_s;
  std::unordered_map<int, std::vector<double>> sigma_t;
  std::unordered_map<int, DealiiMatrix> sigma_t_matrices;
  for (int group = 0; group < n_groups; ++group) {
    sigma_s[group] = MakeMatrix(test_helpers::RandomVector(n_groups*n_groups, 0, max_cross_section_value));
    sigma_t[group] = test_helpers::RandomVector(n_groups, 0, max_cross_section_value);
    DealiiMatrix sigma_t_matrix(n_groups, n_groups);
    for (int i = 0; i < n_groups; ++i)
      sigma_t_matrix(i, i) = sigma_t.at(group).at(i);
    sigma_t_matrices[group] = sigma_t_matrix;
    material_id_mapped_to_spectral_shape_[group] = test_helpers::RandomVector(n_groups, 0, max_spectral_shape_value);
  }

  // set expectations
  ON_CALL(*cross_sections_ptr_, sigma_s()).WillByDefault(Return(sigma_s));
  ON_CALL(*cross_sections_ptr_, sigma_t()).WillByDefault(Return(sigma_t));
  for (int group = 0; group < n_groups; ++group) {
    ON_CALL(*spectral_shape_calculator_obs_ptr_, CalculateSpectralShape(sigma_t_matrices.at(group),
                                                                        sigma_s.at(group)))
        .WillByDefault(Return(material_id_mapped_to_spectral_shape_.at(group)));
  }
  test_calculator_ = std::make_unique<MaterialSpectralShapes>(std::move(spectral_shape_calculator_ptr));
}

auto CalculatorTwoGridMaterialSpectralShapes::MakeMatrix(const std::vector<double> values) const -> DealiiMatrix {
  EXPECT_EQ(values.size(), n_groups * n_groups);
  DealiiMatrix return_matrix(n_groups, n_groups);

  for (int i = 0; i < n_groups; ++i) {
    for (int j = 0; j < n_groups; ++j) {
      return_matrix(i, j) = values.at(i * n_groups + j);
    }
  }

  return return_matrix;
}

TEST_F(CalculatorTwoGridMaterialSpectralShapes, ConstructorAndGetter) {
  ASSERT_NE(test_calculator_->spectral_shape_caluculator_ptr(), nullptr);
  EXPECT_EQ(test_calculator_->spectral_shape_caluculator_ptr(), spectral_shape_calculator_obs_ptr_);
}

TEST_F(CalculatorTwoGridMaterialSpectralShapes, ConstructorThrowsOnNullDependency) {
  EXPECT_ANY_THROW(MaterialSpectralShapes(nullptr));
}

TEST_F(CalculatorTwoGridMaterialSpectralShapes, CalculateSpectralShapes) {
  EXPECT_CALL(*cross_sections_ptr_, sigma_s()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_ptr_, sigma_t()).WillOnce(DoDefault());
  EXPECT_CALL(*spectral_shape_calculator_obs_ptr_, CalculateSpectralShape(_,_))
      .Times(n_groups).WillRepeatedly(DoDefault());

  test_calculator_->CalculateMaterialSpectralShapes(this->cross_sections_ptr_);
  ASSERT_EQ(test_calculator_->material_spectral_shapes().size(), n_groups);
  for (int i = 0; i < n_groups; ++i) {
    EXPECT_EQ(test_calculator_->material_spectral_shapes().at(i), material_id_mapped_to_spectral_shape_.at(i));
  }
}

} // namespace