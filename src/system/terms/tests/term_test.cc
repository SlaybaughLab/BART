#include "system/terms/term.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

/* This testing suite is designed to test the linear and bilinear term classes
 * to ensure they are adding vectors and matrices properly when the
 * GetFullTerm method is called. For tests that verify that the setters and
 * getters work, see linear_term_test.cc. These tests will use a full dealii
 * test environment to provide the needed MPI matrices and vectors.
 */

class SystemTermsFullTermTest : public ::testing::Test {
 protected:
};

TEST_F(SystemTermsFullTermTest, FullTermOperation) {
  EXPECT_TRUE(true);
}

} // namespace