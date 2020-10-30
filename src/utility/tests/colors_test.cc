#include "utility/colors.hpp"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace utility = bart::utility;
namespace test_helpers = bart::test_helpers;

using ::testing::ContainsRegex;

TEST(UtilityColorFreeFunctionTests, ToStringTest) {
  auto test_color = static_cast<utility::Color>(test_helpers::RandomDouble(1, 5));
  std::string regex{"\033\\[3.m"};
  EXPECT_THAT(utility::to_string(test_color),
              ContainsRegex(regex));
}

} // namespace
