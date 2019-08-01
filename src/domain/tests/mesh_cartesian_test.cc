#include "domain/mesh_cartesian.h"

#include <cstdlib>
#include <numeric>
#include <functional>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class DomainMeshCartesianTest : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(DomainMeshCartesianTest, bart::testing::AllDimensions);

TYPED_TEST(DomainMeshCartesianTest, FillTriangulationTest) {
  constexpr int dim = this->dim;

  std::vector<double> spatial_max{btest::RandomVector(dim, 0, 100)};
  std::vector<double> n_cells_double{btest::RandomVector(dim, 1, 20)};
  std::vector<int> n_cells{n_cells_double.begin(), n_cells_double.end()};

  int n_total_cells = std::accumulate(n_cells.begin(), n_cells.end(), 1,
                                      std::multiplies<int>());

  domain::MeshCartesian<dim> test_mesh(spatial_max, n_cells);
  dealii::Triangulation<dim> test_triangulation;

  test_mesh.FillTriangulation(test_triangulation);
  EXPECT_EQ(test_triangulation.n_cells(), n_total_cells);
  EXPECT_FALSE(test_mesh.has_material_mapping());
  for (auto const &cell : test_triangulation.active_cell_iterators()) {
    for (int i = 0; i < dim; ++i) {
      EXPECT_THAT(cell->extent_in_direction(i),
                  ::testing::DoubleNear(spatial_max[i]/n_cells[i], 1e-10));
    }
  }
}

TYPED_TEST(DomainMeshCartesianTest, BadSpatialSize) {
  constexpr int dim = this->dim;

  std::vector<std::vector<double>> spatial_maxes{
      {},
      btest::RandomVector(1, 0, 100),
      btest::RandomVector(2, 0, 100),
      btest::RandomVector(3, 0, 100),
      btest::RandomVector(4, 0, 100)};

  std::vector<std::vector<int>> n_cells{{}, {10}, {10, 20}, {10, 20, 30},
                                        {10, 20, 30, 40}};

  for (int i = -1; i <= 1; i += 2) {
    EXPECT_ANY_THROW({
                       domain::MeshCartesian<dim> test_mesh(spatial_maxes.at(dim + i),
                                                            n_cells.at(dim + i));
                     });
  }
}

} // namespace


class MeshCartesianTest : public ::testing::Test {
 protected:
};

  
class MaterialMapping1DTest : public MeshCartesianTest {
 protected:
  std::vector<double> spatial_max{btest::RandomVector(1, 5, 20)};
  std::vector<int> n_cells{rand() % 20 + 1};
  dealii::Triangulation<1> test_triangulation;
};

TEST_F(MaterialMapping1DTest, 1DMaterialMapping) {
  bart::domain::MeshCartesian<1> test_mesh(spatial_max, n_cells);
  std::string material_mapping{'1'};

  test_mesh.ParseMaterialMap(material_mapping);
  EXPECT_TRUE(test_mesh.has_material_mapping());
  
  std::vector<std::array<double, 1>> test_locations;
  for (int i = 0; i < 5; ++i)
    test_locations.push_back({btest::RandomDouble(0, spatial_max[0])});

  for (auto location : test_locations)
    EXPECT_EQ(test_mesh.GetMaterialID(location), 1) <<
        "Location: " << location[0] << "/" << spatial_max[0];
}

TEST_F(MaterialMapping1DTest, 1DMultiMaterialMapping) {
  bart::domain::MeshCartesian<1> test_mesh(spatial_max, n_cells);
  std::string material_mapping{"1 2"};

  test_mesh.ParseMaterialMap(material_mapping);
  EXPECT_TRUE(test_mesh.has_material_mapping());

  std::vector<std::array<double, 1>> test_locations_1;
  std::vector<std::array<double, 1>> test_locations_2;
  for (int i = 0; i < 5; ++i) {
    test_locations_1.push_back({btest::RandomDouble(0, spatial_max[0]/2)});
    test_locations_2.push_back({btest::RandomDouble(spatial_max[0]/2,
                                                    spatial_max[0])});
  }
  test_locations_1.push_back({0});
  test_locations_1.push_back({spatial_max[0]/2});
  test_locations_2.push_back({spatial_max[0]});

  for (auto location : test_locations_1)
    EXPECT_EQ(test_mesh.GetMaterialID(location), 1) <<
        "Location: " << location[0] << "/" << spatial_max[0];
  for (auto location : test_locations_2)
    EXPECT_EQ(test_mesh.GetMaterialID(location), 2) <<
        "Location: " << location[0] << "/" << spatial_max[0];
}

class MaterialMapping2DTest : public MeshCartesianTest {
 protected:
  std::vector<double> spatial_max{btest::RandomVector(2, 5, 20)};
  std::vector<int> n_cells{rand() % 20 + 1, rand() % 20 + 1};
  dealii::Triangulation<2> test_triangulation;
};

TEST_F(MaterialMapping2DTest, 2DMaterialMapping) {
  bart::domain::MeshCartesian<2> test_mesh(spatial_max, n_cells);
  std::string material_mapping{'1'};

  test_mesh.ParseMaterialMap(material_mapping);
  EXPECT_TRUE(test_mesh.has_material_mapping());

  std::vector<std::array<double, 2>> test_locations;
  for (int i = 0; i < 5; ++i)
    test_locations.push_back({btest::RandomDouble(0, spatial_max[0]),
            btest::RandomDouble(0, spatial_max[1])});

  for (auto location : test_locations)
    EXPECT_EQ(test_mesh.GetMaterialID(location), 1) <<
        "Location: (" << location[0] << ", " << location[1] << ")/("
                      << spatial_max[0] << ", " << spatial_max[1] << ")";
}

TEST_F(MaterialMapping2DTest, 2DMultiMaterialMapping) {
  bart::domain::MeshCartesian<2> test_mesh(spatial_max, n_cells);
  std::string material_mapping{"1 2"};

  test_mesh.ParseMaterialMap(material_mapping);
  EXPECT_TRUE(test_mesh.has_material_mapping());

  std::vector<std::array<double, 2>> test_locations_1;
  std::vector<std::array<double, 2>> test_locations_2;
  for (int i = 0; i < 5; ++i) {
    test_locations_1.push_back(
        {btest::RandomDouble(0, spatial_max[0]/2),          
              btest::RandomDouble(0, spatial_max[1])});
    
    test_locations_2.push_back(
        {btest::RandomDouble(spatial_max[0]/2, spatial_max[0]),
              btest::RandomDouble(0, spatial_max[1])});
  }
  test_locations_1.push_back({0,0});
  test_locations_1.push_back({spatial_max[0]/2, 0});
  test_locations_1.push_back({spatial_max[0]/2, spatial_max[1]/2});
  test_locations_1.push_back({spatial_max[0]/2, spatial_max[1]});
  test_locations_2.push_back({spatial_max[0], 0});
  test_locations_2.push_back({spatial_max[0], spatial_max[1]});

  for (auto location : test_locations_1)
    EXPECT_EQ(test_mesh.GetMaterialID(location), 1) <<
        "Location: (" << location[0] << ", " << location[1] << ")/("
                      << spatial_max[0] << ", " << spatial_max[1] << ")";
  for (auto location : test_locations_2)
    EXPECT_EQ(test_mesh.GetMaterialID(location), 2) <<
        "Location: (" << location[0] << ", " << location[1] << ")/("
                      << spatial_max[0] << ", " << spatial_max[1] << ")";
}

TEST_F(MaterialMapping2DTest, 2DMultiYMaterialMapping) {
  bart::domain::MeshCartesian<2> test_mesh(spatial_max, n_cells);
  std::string material_mapping{"2 1\n1 2"};

  test_mesh.ParseMaterialMap(material_mapping);
  EXPECT_TRUE(test_mesh.has_material_mapping());

  std::vector<std::array<double, 2>> test_locations_1;
  std::vector<std::array<double, 2>> test_locations_2;
  for (int i = 0; i < 5; ++i) {
    test_locations_1.push_back(
        {btest::RandomDouble(0, spatial_max[0]/2),          
              btest::RandomDouble(0, spatial_max[1]/2)});
    
    test_locations_1.push_back(
        {btest::RandomDouble(spatial_max[0]/2, spatial_max[0]),
              btest::RandomDouble(spatial_max[1]/2, spatial_max[1])});
    
    test_locations_2.push_back(
        {btest::RandomDouble(spatial_max[0]/2, spatial_max[0]),
              btest::RandomDouble(0, spatial_max[1]/2)});
    
    test_locations_2.push_back(
        {btest::RandomDouble(0, spatial_max[0]/2),
              btest::RandomDouble(spatial_max[1]/2, spatial_max[1])});
  }
  test_locations_1.push_back({0,0});
  test_locations_1.push_back({spatial_max[0]/2, 0});
  test_locations_1.push_back({spatial_max[0]/2, spatial_max[1]/2});
  test_locations_2.push_back({spatial_max[0]/2, spatial_max[1]});
  test_locations_2.push_back({spatial_max[0], 0});
  test_locations_1.push_back({spatial_max[0], spatial_max[1]});
  test_locations_2.push_back({spatial_max[0]/2, spatial_max[1]});

  for (auto location : test_locations_1)
    EXPECT_EQ(test_mesh.GetMaterialID(location), 1) <<
        "Location: (" << location[0] << ", " << location[1] << ")/("
                      << spatial_max[0] << ", " << spatial_max[1] << ")";
  for (auto location : test_locations_2)
    EXPECT_EQ(test_mesh.GetMaterialID(location), 2) <<
        "Location: (" << location[0] << ", " << location[1] << ")/("
                      << spatial_max[0] << ", " << spatial_max[1] << ")";
}
