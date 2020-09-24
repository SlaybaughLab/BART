#include "instrumentation/converter/pair_incrementer.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace converter = bart::instrumentation::converter;
namespace test_helpers = bart::test_helpers;

class InstrumentationConverterIncrementerTest : public ::testing::Test {
 public:
};

TEST_F(InstrumentationConverterIncrementerTest, Increment) {
  using Incrementer = converter::PairIncrementer<double>;

  double first_input{test_helpers::RandomDouble(-100, 100)};
  double second_input{test_helpers::RandomDouble(-100, 100)};

  Incrementer test_incrementer;

  auto first_pair = test_incrementer.Convert(first_input);
  EXPECT_THAT(first_pair, ::testing::Pair(0, first_input));
  auto second_pair = test_incrementer.Convert(second_input);
  EXPECT_THAT(second_pair, ::testing::Pair(1, second_input));
}

} // namespace
