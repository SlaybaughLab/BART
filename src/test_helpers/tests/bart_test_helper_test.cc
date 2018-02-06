#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <sys/stat.h>
#include <exception>

#include "../bart_test_helper.h"

using ::testing::MatchesRegex;

class BartTestHelperTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  std::string gold_files_directory = "test_data/";
};

TEST_F(BartTestHelperTest, InitializationNoReport) {
  btest::BartTestHelper test_helper;
  ASSERT_EQ(test_helper.GetReportDirectory(), "");
}

TEST_F(BartTestHelperTest, InitalizeReport) {
  btest::BartTestHelper test_helper(true, gold_files_directory);
  std::string report_directory = test_helper.GetReportDirectory();
  // Verify report directory name
  EXPECT_THAT(report_directory,
              MatchesRegex(gold_files_directory + "........_...._fail"));
  // Verify report directory existence
  struct stat sb;
  ASSERT_TRUE(stat(report_directory.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode));
  rmdir(report_directory.c_str());
}

TEST_F(BartTestHelperTest, InitalizeBadDirectory) {
  std::string bad_directory = "testing_data/";
  ASSERT_THROW(btest::BartTestHelper test_helper(true, bad_directory),
               std::runtime_error);
}
  
