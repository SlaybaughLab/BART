#include "formulation/scalar/cfem_diffusion.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class FormulationCFEMDiffusionTest : public ::testing::Test {};

TEST_F(FormulationCFEMDiffusionTest, ConstructorTest) {
  formulation::scalar::CFEM_Diffusion test_diffusion();
}

} // namespace