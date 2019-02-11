#include "../cartesian_mesh.h"

#include <cstdlib>
#include <numeric>
#include <functional>
#include <vector>

#include <deal.II/grid/tria.h>
#include <gtest/gtest.h>

#include "../../test_helpers/test_helper_functions.h"
#include "../../test_helpers/gmock_wrapper.h"

class CartesianMeshTest : public ::testing::Test {
 protected:
  template <int dim> void FillTriangulationTest();
};

template <int dim>
void CartesianMeshTest::FillTriangulationTest() {
  std::vector<double> spatial_max{btest::RandomVector(dim, 0, 100)};
  std::vector<double> n_cells_double{btest::RandomVector(dim, 0, 20)};
  std::vector<int> n_cells{n_cells_double.begin(), n_cells_double.end()};

  int n_total_cells = std::accumulate(n_cells.begin(), n_cells.end(), 1,
                                      std::multiplies<int>());

  bart::domain::CartesianMesh<dim> test_mesh(spatial_max, n_cells);
  dealii::Triangulation<dim> test_triangulation;

  test_mesh.FillTriangulation(test_triangulation);
  EXPECT_EQ(test_triangulation.n_cells(), n_total_cells);

  for (auto const &cell : test_triangulation.active_cell_iterators()) {
    for (int i = 0; i < dim; ++i) {
      EXPECT_THAT(cell->extent_in_direction(i),
                  ::testing::DoubleNear(spatial_max[i]/n_cells[i], 1e-14));
    }
  }
}                                     

TEST_F(CartesianMeshTest, Triangulation1D) {
    FillTriangulationTest<1>();
}

TEST_F(CartesianMeshTest, Triangulation2D) {
    FillTriangulationTest<2>();
}

TEST_F(CartesianMeshTest, Triangulation3D) {
    FillTriangulationTest<3>();
}

TEST_F(CartesianMeshTest, BadSpatialSize) {

  std::vector<std::vector<double>> domain_sizes {{10.0}, {10.0, 20.0},
                                                         {10.0, 20.0, 30.0}};
  std::vector<std::vector<int>> n_cells {{10}, {10, 20}, {10, 20, 30}};
  
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      // 1D Case
      if ((i != 0) || (j != 0)) {
        EXPECT_ANY_THROW({
            bart::domain::CartesianMesh<1> test_mesh(domain_sizes[i], n_cells[j]);
          });
      } else {
        EXPECT_NO_THROW({
            bart::domain::CartesianMesh<1> test_mesh(domain_sizes[i], n_cells[j]);
          });
      }
      // 2D Case
      if ((i != 1) || (j != 1)) {
        EXPECT_ANY_THROW({
            bart::domain::CartesianMesh<2> test_mesh(domain_sizes[i], n_cells[j]);
          });
      } else {
        EXPECT_NO_THROW({
            bart::domain::CartesianMesh<2> test_mesh(domain_sizes[i], n_cells[j]);
          });
      }
      // 3D Case
      if ((i != 2) || (j != 2)) {
        EXPECT_ANY_THROW({
            bart::domain::CartesianMesh<3> test_mesh(domain_sizes[i], n_cells[j]);
          });
      } else {
        EXPECT_NO_THROW({
            bart::domain::CartesianMesh<3> test_mesh(domain_sizes[i], n_cells[j]);
          });
      }
    }
  }
}

class MaterialMappingTest : public CartesianMeshTest {
 protected:
  template <int dim> std::array<double, dim> RandomArray(int min, int max);
};

template <int dim>
std::array<double, dim> MaterialMappingTest::RandomArray(int min, int max) {
  auto vector = btest::RandomVector(dim, min, max);
  std::array<double, dim> return_array;
  std::copy(vector.begin(), vector.end(), return_array.begin());
  return return_array;
}
  
class MaterialMapping1DTest : public MaterialMappingTest {
 protected:
  std::vector<double> spatial_max{btest::RandomVector(1, 0, 20)};
  std::vector<int> n_cells{rand() % 20 + 1};
  dealii::Triangulation<1> test_triangulation;
};

TEST_F(MaterialMapping1DTest, 1DMaterialMapping) {
  bart::domain::CartesianMesh<1> test_mesh(spatial_max, n_cells);
  std::string material_mapping{'1'};

  test_mesh.SetMaterialIDs(test_triangulation, material_mapping);

  std::vector<std::array<double, 1>> test_locations
  {
    RandomArray<1>(0, spatial_max[0]),
        RandomArray<1>(0, spatial_max[0]),
        RandomArray<1>(0, spatial_max[0])};

  for (auto location : test_locations)
    ASSERT_EQ(test_mesh.GetMaterialID(location), 1);
}

template std::array<double, 1> MaterialMappingTest::RandomArray<1>(int, int);
template std::array<double, 2> MaterialMappingTest::RandomArray<2>(int, int);
template std::array<double, 3> MaterialMappingTest::RandomArray<3>(int, int);

template void CartesianMeshTest::FillTriangulationTest<1>();
template void CartesianMeshTest::FillTriangulationTest<2>();
template void CartesianMeshTest::FillTriangulationTest<3>();
