#include "../domain.h"

#include <memory>

#include <gtest/gtest.h>

#include "../../test_helpers/gmock_wrapper.h"
#include "../../domain/tests/mesh_mock.h"
#include "../../domain/tests/finite_element_mock.h"

class DomainTest : public ::testing::Test {
 protected:

  std::unique_ptr<bart::domain::MeshI<2>> mesh_ptr;
  std::unique_ptr<bart::domain::FiniteElementI<2>> fe_ptr;
  
  void SetUp() override;
};

void DomainTest::SetUp() {
  mesh_ptr = std::make_unique<bart::domain::MeshMock<2>>();
  fe_ptr = std::make_unique<bart::domain::FiniteElementMock<2>>();
}

TEST_F(DomainTest, Constructor) {
  bart::problem::Domain<2> test_domain(mesh_ptr, fe_ptr);
  EXPECT_EQ(mesh_ptr, nullptr);
  EXPECT_EQ(fe_ptr, nullptr);
}

