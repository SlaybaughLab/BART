#include "acceleration/two_grid/spectral_shape/domain_spectral_shapes.hpp"

#include "domain/tests/domain_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_helper_functions.h"
#include "acceleration/two_grid/spectral_shape/tests/domain_spectral_shapes_mock.hpp"

namespace  {

using namespace bart;
using ::testing::Return, ::testing::ContainerEq;

template <typename DimensionWrapper>
class AccelerationTwoGridDomainSpectralShapesTest : public bart::testing::DealiiTestDomain<DimensionWrapper::value>,
                                                    public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  AccelerationTwoGridDomainSpectralShapesTest() : bart::testing::DealiiTestDomain<dim>(1.0, 2) {};
  using DomainMock = domain::DomainMock<dim>;
  using Vector = dealii::Vector<double>;
  using DomainSpectralShapeMap = std::unordered_map<int, Vector>; // Mapping of group to spectral shape value by dof
  using MaterialSpectralShapeMap = std::unordered_map<int, std::vector<double>>; // Mapping of material ID to spectral shape value by group

  DomainSpectralShapeMap expected_domain_spectral_shape_map_;
  MaterialSpectralShapeMap material_spectral_shape_map_;

  acceleration::two_grid::spectral_shape::DomainSpectralShapes<dim> test_calculator;
  std::shared_ptr<DomainMock> domain_mock_ptr_{ std::make_shared<DomainMock>() };

  // Test parameters
  const int n_groups{ test_helpers::RandomInt(1, 10) };
  const int n_materials{ test_helpers::RandomInt(1, 10) };

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto AccelerationTwoGridDomainSpectralShapesTest<DimensionWrapper>::SetUp() -> void {
  this->SetUpDealii();

  // Set up the two mappings
  for (int material_id = 0; material_id < n_materials; ++material_id)
    material_spectral_shape_map_[material_id] = test_helpers::RandomVector(n_groups, 0, 10);
  for (int group = 0; group < n_groups; ++group)
    expected_domain_spectral_shape_map_[group] = Vector(this->dof_handler_.n_dofs());

  // Assign random material values and then fill in expected domain spectral shape based on cell DOFS
  for (auto& cell : this->cells_) {
    const int material_id { test_helpers::RandomInt(0, n_materials) };
    cell->set_material_id(material_id);
    std::vector<unsigned int> global_dofs(this->fe_.dofs_per_cell);
    cell->get_dof_indices(global_dofs);
    for (int group = 0; group < n_groups; ++group) {
      for (int index : global_dofs) {
        expected_domain_spectral_shape_map_[group][index] += material_spectral_shape_map_[material_id].at(group);
      }
    }
  }
}


TYPED_TEST_SUITE(AccelerationTwoGridDomainSpectralShapesTest, bart::testing::AllDimensions);

TYPED_TEST(AccelerationTwoGridDomainSpectralShapesTest, CalculateDomainSpectralShape) {
  dealii::Vector<double> cell_vector(this->fe_.dofs_per_cell);
  EXPECT_CALL(*this->domain_mock_ptr_, Cells()).WillOnce(Return(this->cells_));
  EXPECT_CALL(*this->domain_mock_ptr_, total_degrees_of_freedom()).WillOnce(Return(this->dof_handler_.n_dofs()));
  EXPECT_CALL(*this->domain_mock_ptr_, GetCellVector()).WillOnce(Return(cell_vector));
  const auto domain_spectral_shape_map = this->test_calculator.CalculateDomainSpectralShapes(
      this->material_spectral_shape_map_, *this->domain_mock_ptr_);
  ASSERT_EQ(domain_spectral_shape_map.size(), this->n_groups);
  EXPECT_THAT(domain_spectral_shape_map, ContainerEq(this->expected_domain_spectral_shape_map_));
}

} // namespace
