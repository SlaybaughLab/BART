#include "../cartesian_mesh.h"

#include <array>
#include <deal.II/grid/tria.h>
#include <gtest/gtest.h>

#include "../../test_helpers/gmock_wrapper.h"

class CartesianMeshTest : public ::testing::Test {
 protected:
};

TEST_F(CartesianMeshTest, Triangulation1D) {
  bart::domain::CartesianMesh<1> test_mesh;
  dealii::Triangulation<1> test_triangulation;
  std::array<double, 1> x_max{10.5};
  std::array<int, 1> n_cells{10};
  test_mesh.FillTriangulation(test_triangulation, x_max, n_cells);
  
  EXPECT_EQ(test_triangulation.n_lines(), 10);
  EXPECT_EQ(test_triangulation.n_cells(), 10);
  EXPECT_EQ(test_triangulation.n_used_vertices(), 11);

  for (auto const &cell : test_triangulation.active_cell_iterators()) {
    EXPECT_THAT(cell->extent_in_direction(0),
                ::testing::DoubleNear(x_max[0]/n_cells[0], 1e-14));
  }
}

TEST_F(CartesianMeshTest, Triangulation2D) {
  bart::domain::CartesianMesh<2> test_mesh;
  dealii::Triangulation<2> test_triangulation;
  std::array<double, 2> max_dim{10.5, 12.0};
  std::array<int, 2> n_cells{12, 20};
  test_mesh.FillTriangulation(test_triangulation, max_dim, n_cells);
  
  EXPECT_EQ(test_triangulation.n_cells(), 240);

  for (auto const &cell : test_triangulation.active_cell_iterators()) {
    for (int i = 0; i < 2; ++i) {
      EXPECT_THAT(cell->extent_in_direction(i),
                  ::testing::DoubleNear(max_dim[i]/n_cells[i], 1e-14));
    }
  }
}

TEST_F(CartesianMeshTest, Triangulation3D) {
  bart::domain::CartesianMesh<3> test_mesh;
  dealii::Triangulation<3> test_triangulation;
  std::array<double, 3> max_dim{10.5, 12.0, 6.0};
  std::array<int, 3> n_cells{12, 20, 10};
  test_mesh.FillTriangulation(test_triangulation, max_dim, n_cells);
  
  EXPECT_EQ(test_triangulation.n_cells(), 12*20*10);

  for (auto const &cell : test_triangulation.active_cell_iterators()) {
    for (int i = 0; i < 3; ++i) {
      EXPECT_THAT(cell->extent_in_direction(i),
                  ::testing::DoubleNear(max_dim[i]/n_cells[i], 1e-14));
    }
  }
}
