#include "../domain.h"

#include <gtest/gtest.h>

#include "../../test_helpers/gmock_wrapper.h"
#include "../../domain/tests/mesh_mock.h"
#include "../../domain/tests/finite_element_mock.h"

class DomainTest : public ::testing::Test {
 protected:
  bart::domain::MeshMock<2> mock_mesh;
  bart::domain::FiniteElementMock<2> mock_finite_element;
};

TEST_F(DomainTest, ProvideTest) {
  
}
