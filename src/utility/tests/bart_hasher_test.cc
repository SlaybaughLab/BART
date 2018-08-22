#include "../bart_hasher.h"

#include <gtest/gtest.h>
#include <string>

class BartHasherTest : public ::testing::Test {
 protected:
  butil::Hasher test_hasher;
};

TEST_F(BartHasherTest, ArrayTest) {
  std::array<int, 2> test_array{1,2};
  std::size_t result(test_hasher(test_array));
  std::string hash{"3370697991563800380"};
  EXPECT_EQ(std::to_string(result), hash);
}

TEST_F(BartHasherTest, VectorTest) {
  std::vector<int> test_vector{1,2};
  std::size_t result(test_hasher(test_vector));
  std::string hash{"3370697991563800380"};
  EXPECT_EQ(std::to_string(result), hash);
}
