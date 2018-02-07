#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <sys/stat.h>
#include <exception>
#include <fstream>

#include "../bart_test_helper.h"

using ::testing::MatchesRegex;

class TestHelperIntTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  std::string gold_files_directory = "test_data/";
};

TEST_F(TestHelperIntTest, IntegrationTestGoodNoReport) {
  btest::BartTestHelper test_helper(false, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper";
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "1234";
  actual_stream.close();
  EXPECT_TRUE(test_helper.GoldTest(filename));
  // Check deleted
  struct stat sb;
  EXPECT_FALSE(stat(filename.c_str(), &sb) == 0);
}

TEST_F(TestHelperIntTest, IntegrationTestBadNoReport) {
  btest::BartTestHelper test_helper(false, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper";
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "1235";
  actual_stream.close();
  EXPECT_FALSE(test_helper.GoldTest(filename));
  // Check deleted
  struct stat sb;
  EXPECT_FALSE(stat(filename.c_str(), &sb) == 0);
}

TEST_F(TestHelperIntTest, IntegrationTestNoActual) {
  btest::BartTestHelper test_helper(false, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper";
  EXPECT_FALSE(test_helper.GoldTest(filename));
}

TEST_F(TestHelperIntTest, IntegrationTestBadGoldFileNoReport) {
  std::string bad_gold_dir = "/testing_files";
  btest::BartTestHelper test_helper(false, bad_gold_dir);
  // Make actual file
  std::string filename = "bart_test_helper";
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "1235";
  actual_stream.close();
  EXPECT_FALSE(test_helper.GoldTest(filename));
  // Check deleted
  struct stat sb;
  EXPECT_FALSE(stat(filename.c_str(), &sb) == 0);
}

TEST_F(TestHelperIntTest, IntegrationTestBadGoldReport) {
  btest::BartTestHelper test_helper(true, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper_2";
  std::string report_directory = test_helper.GetReportDirectory();
  
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "1235";
  actual_stream.close();
  EXPECT_FALSE(test_helper.GoldTest(filename));
  // Check deleted
  struct stat sb;
  EXPECT_FALSE(stat(filename.c_str(), &sb) == 0);
    std::string new_name(report_directory + "/" + filename);
  EXPECT_TRUE(
      stat(new_name.c_str(), &sb) == 0 &&
      S_ISREG(sb.st_mode));
  EXPECT_FALSE(
      stat((new_name + ".diff").c_str(), &sb) == 0 &&
      S_ISREG(sb.st_mode));
  remove(new_name.c_str());
  rmdir(report_directory.c_str());
}

TEST_F(TestHelperIntTest, IntegrationTestBadReport) {
  btest::BartTestHelper test_helper(true, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper";
  std::string report_directory = test_helper.GetReportDirectory();
  
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "12\n34";
  actual_stream.close();
  EXPECT_FALSE(test_helper.GoldTest(filename));
  // Check deleted
  struct stat sb;
  std::string new_name(report_directory + "/" + filename);
  std::string diff_name = new_name + ".diff";
  //Files in correct places
  EXPECT_FALSE(stat(filename.c_str(), &sb) == 0);
  EXPECT_TRUE(
      stat(new_name.c_str(), &sb) == 0 &&
      S_ISREG(sb.st_mode));
  EXPECT_TRUE(
      stat(diff_name.c_str(), &sb) == 0 &&
      S_ISREG(sb.st_mode));

  // Cleanup
  remove(new_name.c_str());
  remove(diff_name.c_str());
  rmdir(report_directory.c_str());
}
