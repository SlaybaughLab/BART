#include "data/system/term.h"

#include <memory>

#include "data/system/system_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using VariableLinearTerms = data::system::VariableLinearTerms;

class SystemLinearTermTest : public ::testing::Test {
 protected:

  data::system::Term<data::system::MPILinearTermPair> test_linear_term;
};

TEST_F(SystemLinearTermTest, Dummy) {
  EXPECT_TRUE(true);
}


} // namespace

