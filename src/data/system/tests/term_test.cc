#include "data/system/term.h"

#include "data/system/system_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename TermPair>
class SystemTermTest : public ::testing::Test {
 protected:

  data::system::Term<TermPair> test_term_;
};

using TermTypes = ::testing::Types<
    data::system::MPILinearTermPair,
    data::system::MPIBilinearTermPair
>;

TYPED_TEST_CASE(SystemTermTest, TermTypes);

TYPED_TEST(SystemTermTest, Dummy) {
  EXPECT_TRUE(true);
}

} // namespace

