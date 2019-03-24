#include "convergence/parameters/single_parameter_checker.h"

#include "convergence/tests/single_checker_test.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

class ConvergenceSingleParameterCheckerTest :
 public bart::convergence::testing::SingleCheckerTest<double>{

 protected:
  bart::convergence::parameters::SingleParameterChecker test_checker{1e-6};

};

TEST_F(ConvergenceSingleParameterCheckerTest, BaseClassTests) {
  TestBaseMethods(&test_checker);
}


} // namespace

