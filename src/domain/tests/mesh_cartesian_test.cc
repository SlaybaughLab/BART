#include "domain/mesh_cartesian.h"

#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <functional>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>
#include <problem/parameter_types.h>

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

TYPED_TEST(DomainMeshCartesianTest, FillBoundaryIDTest) {
  constexpr int dim = this->dim;
  using Boundary = problem::Boundary;

  std::vector<double> spatial_max(dim, 10.0);
  std::vector<int> n_cells(dim, 2);

  const double zero_tol = 1.0e-14;

  domain::MeshCartesian<dim> test_mesh(spatial_max, n_cells);
  dealii::Triangulation<dim> test_triangulation;

  test_mesh.FillTriangulation(test_triangulation);
  test_mesh.FillBoundaryID(test_triangulation);

  int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;

  for (auto const &cell : test_triangulation.active_cell_iterators()) {
    if (cell->is_locally_owned()) {
      for (int face_id = 0; face_id < faces_per_cell; ++face_id) {
        auto face = cell->face(face_id);
        dealii::Point<dim> face_center = face->center();
        switch (dim) {
          case 3: {
            if (std::fabs(face_center[2]) < zero_tol) {
              EXPECT_EQ(static_cast<Boundary>(face->boundary_id()),Boundary::kZMin);
              break;
            } else if (std::fabs(face_center[2] - spatial_max.at(2)) < zero_tol) {
              EXPECT_EQ(static_cast<Boundary>(face->boundary_id()), Boundary::kZMax);
              break;
            }
            [[fallthrough]];
          }
          case 2: {
            if (std::fabs(face_center[1]) < zero_tol) {
              EXPECT_EQ(static_cast<Boundary>(face->boundary_id()), Boundary::kYMin);
              break;
            } else if (std::fabs(face_center[1] - spatial_max.at(1)) < zero_tol) {
              EXPECT_EQ(static_cast<Boundary>(face->boundary_id()), Boundary::kYMax);
              break;
            }
            [[fallthrough]];
          }
            // Fall through to check x-direction
          case 1: {
            if (std::fabs(face_center[0]) < zero_tol) {
              EXPECT_EQ(static_cast<Boundary>(face->boundary_id()), Boundary::kXMin);
              break;
            } else if (std::fabs(face_center[0] - spatial_max.at(0)) < zero_tol) {
              EXPECT_EQ(static_cast<Boundary>(face->boundary_id()), Boundary::kXMax);
              break;
            }
            break;
          }
          default: {
            AssertThrow(false,
                        dealii::ExcMessage("Unsupported number of dimensions in FillBoundaryID"));
          }
        }
      }
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
  std::array<int, 2> i_values{-1, 1};

  for (const auto& i : i_values) {
    EXPECT_ANY_THROW({
                       domain::MeshCartesian<dim> test_mesh(spatial_maxes.at(dim + i),
                                                            n_cells.at(dim + i));
                     });
  }
}

TYPED_TEST(DomainMeshCartesianTest, SingleMaterialMapping) {
  constexpr int dim = this->dim;
  std::vector<double> spatial_max{btest::RandomVector(dim, 5, 20)};
  std::vector<double> n_cells_double{btest::RandomVector(dim, 5, 20)};
  std::vector<int> n_cells{n_cells_double.cbegin(),
                           n_cells_double.cend()};

  std::string material_mapping{'1'};
  domain::MeshCartesian<dim> test_mesh(spatial_max, n_cells, material_mapping);

  EXPECT_TRUE(test_mesh.has_material_mapping());

  // Random inner locations
  std::array<std::array<double, dim>, 10> test_locations;

  for (auto& location : test_locations) {
    for (int i = 0; i < dim; ++i) {
      location.at(i) = btest::RandomDouble(0, spatial_max.at(i));
    }

    EXPECT_EQ(test_mesh.GetMaterialID(location), 1);
  }

  // Edges and corners
  std::array<double, dim> test_location;
  dealii::Point<dim> test_point;

  std::array<double, 3> x_locations{0, spatial_max.at(0)/2, spatial_max.at(0)};
  std::vector<double> y_locations{};
  std::vector<double> z_locations{};

  if (dim > 1) {
    y_locations = {0, spatial_max.at(1)/2, spatial_max.at(1)};
    if (dim > 2) {
      z_locations = {0, spatial_max.at(2)/2, spatial_max.at(2)};
    }
  }

  for (const int x : x_locations) {
    test_location.at(0) = x;
    test_point[0] = x;
    for (const int y : y_locations) {
      test_location.at(1) = y;
      test_point[1] = y;
      for (const int z : z_locations) {
        test_location.at(2) = z;
        test_point[2] = z;
      }
    }
    EXPECT_EQ(test_mesh.GetMaterialID(test_location), 1);
    EXPECT_EQ(test_mesh.GetMaterialID(test_point), 1);
  }
}

class DomainMeshCartesianMappingTest : public ::testing::Test {};

TEST_F(DomainMeshCartesianMappingTest, MultipleMaterialMapping1D) {
  std::vector<double> spatial_max{btest::RandomVector(1, 5, 20)};
  std::vector<double> n_cells_double{btest::RandomVector(1, 5, 20)};
  std::vector<int> n_cells{n_cells_double.cbegin(), n_cells_double.cend()};

  domain::MeshCartesian<1> test_mesh(spatial_max, n_cells);

  std::string material_mapping{"1 2"};

  test_mesh.ParseMaterialMap(material_mapping);
  EXPECT_TRUE(test_mesh.has_material_mapping());

  std::array<std::array<double, 1>, 5> test_locations;

  for (auto& location : test_locations) {
    location.at(0) = btest::RandomDouble(0, spatial_max.at(0)/2);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 1);
    location.at(0) = btest::RandomDouble(spatial_max.at(0)/2, spatial_max.at(0));
    EXPECT_EQ(test_mesh.GetMaterialID(location), 2);
  }

  // Edge cases
  std::array<double, 1> origin{0};
  std::array<double, 1> midpoint{spatial_max.at(0)/2};
  std::array<double, 1> endpoint{spatial_max.at(0)};

  EXPECT_EQ(test_mesh.GetMaterialID(origin), 1);
  EXPECT_EQ(test_mesh.GetMaterialID(midpoint), 1);
  EXPECT_EQ(test_mesh.GetMaterialID(endpoint), 2);
}

TEST_F(DomainMeshCartesianMappingTest, MultipleMaterialMapping2D) {
  std::vector<double> spatial_max{btest::RandomVector(2, 5, 20)};
  std::vector<double> n_cells_double{btest::RandomVector(2, 5, 20)};
  std::vector<int> n_cells{n_cells_double.cbegin(), n_cells_double.cend()};

  domain::MeshCartesian<2> test_mesh(spatial_max, n_cells);

  std::string material_mapping{"1 2\n3 4"};
  test_mesh.ParseMaterialMap(material_mapping);

  double x_max = spatial_max.at(0), x_mid = spatial_max.at(0)/2;
  double y_max = spatial_max.at(1), y_mid = spatial_max.at(1)/2;

  EXPECT_TRUE(test_mesh.has_material_mapping());
  // Inner locations
  std::array<std::array<double, 2>, 5> test_locations;
  for (auto& location : test_locations) {
    location.at(0) = btest::RandomDouble(0, x_mid);
    location.at(1) = btest::RandomDouble(0, y_mid);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 3);

    location.at(0) = btest::RandomDouble(x_mid, x_max);
    location.at(1) = btest::RandomDouble(0, y_mid);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 4);

    location.at(0) = btest::RandomDouble(0, x_mid);
    location.at(1) = btest::RandomDouble(y_mid, y_max);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 1);

    location.at(0) = btest::RandomDouble(x_mid, x_max);
    location.at(1) = btest::RandomDouble(y_mid, y_max);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 2);
  }
  // Edges and corners
  std::map<std::array<double, 2>, int> locations_and_material_ids{
      {{0,0}, 3}, {{x_mid, 0}, 3}, {{x_max, 0}, 4},
      {{0, y_mid}, 3}, {{x_mid, y_mid}, 3}, {{x_max, y_mid}, 4},
      {{0, y_max}, 1}, {{x_mid, y_max}, 1}, {{x_max, y_max}, 2}
  };

  for (auto& location_and_material_id : locations_and_material_ids) {
    auto& [location, id] = location_and_material_id;
    EXPECT_EQ(test_mesh.GetMaterialID(location), id);
  }

}

TEST_F(DomainMeshCartesianMappingTest, MultipleMaterialMapping3D) {
  std::vector<double> spatial_max{btest::RandomVector(3, 5, 20)};
  std::vector<double> n_cells_double{btest::RandomVector(3, 5, 20)};
  std::vector<int> n_cells{n_cells_double.cbegin(), n_cells_double.cend()};

  domain::MeshCartesian<3> test_mesh(spatial_max, n_cells);

  std::string material_mapping{"1 2\n3 4\n\n5 6\n7 8"};
  test_mesh.ParseMaterialMap(material_mapping);

  double x_max = spatial_max.at(0), x_mid = spatial_max.at(0)/2;
  double y_max = spatial_max.at(1), y_mid = spatial_max.at(1)/2;
  double z_max = spatial_max.at(2), z_mid = spatial_max.at(2)/2;

  EXPECT_TRUE(test_mesh.has_material_mapping());
  // Inner locations
  std::array<std::array<double, 3>, 5> test_locations;
  for (auto& location : test_locations) {

    location.at(0) = btest::RandomDouble(0, x_mid);
    location.at(1) = btest::RandomDouble(0, y_mid);
    location.at(2) = btest::RandomDouble(0, z_mid);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 7);

    location.at(0) = btest::RandomDouble(x_mid, x_max);
    location.at(1) = btest::RandomDouble(0, y_mid);
    location.at(2) = btest::RandomDouble(0, z_mid);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 8);

    location.at(0) = btest::RandomDouble(0, x_mid);
    location.at(1) = btest::RandomDouble(y_mid, y_max);
    location.at(2) = btest::RandomDouble(0, z_mid);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 5);

    location.at(0) = btest::RandomDouble(0, x_mid);
    location.at(1) = btest::RandomDouble(0, y_mid);
    location.at(2) = btest::RandomDouble(z_mid, z_max);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 3);

    location.at(0) = btest::RandomDouble(x_mid, x_max);
    location.at(1) = btest::RandomDouble(y_mid, y_max);
    location.at(2) = btest::RandomDouble(0, z_mid);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 6);

    location.at(0) = btest::RandomDouble(x_mid, x_max);
    location.at(1) = btest::RandomDouble(y_mid, y_max);
    location.at(2) = btest::RandomDouble(z_mid, z_max);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 2);

    location.at(0) = btest::RandomDouble(x_mid, x_max);
    location.at(1) = btest::RandomDouble(0, y_mid);
    location.at(2) = btest::RandomDouble(z_mid, z_max);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 4);

    location.at(0) = btest::RandomDouble(0, x_mid);
    location.at(1) = btest::RandomDouble(y_mid, y_max);
    location.at(2) = btest::RandomDouble(z_mid, z_max);
    EXPECT_EQ(test_mesh.GetMaterialID(location), 1);
  }
}


} // namespace
