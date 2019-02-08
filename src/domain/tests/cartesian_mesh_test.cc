#include "../cartesian_mesh.h"

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
  bart::domain::CartesianMesh<dim> test_mesh;
  dealii::Triangulation<dim> test_triangulation;

  std::vector<double> spatial_max{btest::RandomVector(dim, 0, 100)};
  std::vector<double> n_cells_double{btest::RandomVector(dim, 0, 20)};
  std::vector<int> n_cells{n_cells_double.begin(), n_cells_double.end()};

  int n_total_cells = std::accumulate(n_cells.begin(), n_cells.end(), 1,
                                      std::multiplies<int>());

  test_mesh.FillTriangulation(test_triangulation, spatial_max, n_cells);
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

template void CartesianMeshTest::FillTriangulationTest<1>();
template void CartesianMeshTest::FillTriangulationTest<2>();
template void CartesianMeshTest::FillTriangulationTest<3>();
