#include "../self_adjoint_angular_flux.h"

#include <sstream>

#include <gtest/gtest.h>

#include "../../test_helpers/gmock_wrapper.h"
#include "../../common/problem_definition.h"

class SAAFIntegrationTest : public ::testing::Test {
 protected:
  void SetUp() override;
  dealii::ParameterHandler prm_;
};

void SAAFIntegrationTest::SetUp() {
  bparams::DeclareParameters(prm_);
  std::ostringstream parameter_oss;
  parameter_oss << "set transport model          = saaf\n"
                << "set angular quadrature name  = lsgc\n"
                << "set angular quadrature order = 2\n"
                << "set ho linear solver name    = direct\n"
                << "set ho preconditioner name   = bssor\n";
  const char *parse_string = parameter_oss.str().c_str();
  prm_.parse_input_from_string(parse_string);
}
