#include "domain/mesh/factory.hpp"
#include "domain/mesh/mesh_cartesian.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace mesh = bart::domain::mesh;
namespace test_helpers = bart::test_helpers;
using MeshName = mesh::MeshName;

template <typename DimensionWrapper>
class DomainMeshFactoryTest : public ::testing::Test {
 public:
  constexpr static int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(DomainMeshFactoryTest, bart::testing::AllDimensions);

TYPED_TEST(DomainMeshFactoryTest, MeshInstantiation) {
  constexpr int dim = this->dim;

  std::vector<double> spatial_max{test_helpers::RandomVector(dim, 0, 100)};
  std::vector<double> n_cells_double{test_helpers::RandomVector(dim, 1, 20)};
  std::vector<int> n_cells{n_cells_double.begin(), n_cells_double.end()};
  std::string material_mapping{ "1" };

  auto mesh_ptr = mesh::MeshIFactory<dim, const std::vector<double>, const std::vector<int>, const std::string>::get()
      .GetConstructor(MeshName::kCartesian)(spatial_max, n_cells, material_mapping);
  using ExpectedType = mesh::MeshCartesian<this->dim>;
  ASSERT_NE(mesh_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<ExpectedType*>(mesh_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
}

} // namespace
