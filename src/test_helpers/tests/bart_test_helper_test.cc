#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/logstream.h>

#include <sys/stat.h>
#include <exception>
#include <fstream>

#include "../bart_test_helper.h"

using ::testing::MatchesRegex;

class BartTestHelperTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  std::string gold_files_directory = "test_data/";
  btest::BartTestHelper test_helper;
  virtual void TearDown() {
    std::string report_directory = test_helper.GetReportDirectory();
    if (report_directory != "")
      rmdir(test_helper.GetReportDirectory().c_str());
  }
};

TEST_F(BartTestHelperTest, InitializationNoReport) {
  ASSERT_EQ(test_helper.GetReportDirectory(), "");
}

TEST_F(BartTestHelperTest, ConstructorReport) {
  btest::BartTestHelper test_helper2(true, gold_files_directory);
  std::string report_directory = test_helper2.GetReportDirectory();
  // Verify report directory name
  EXPECT_THAT(report_directory,
              MatchesRegex(gold_files_directory + "........_...._.._report"));
  // Verify report directory existence
  struct stat sb;
  ASSERT_TRUE(stat(report_directory.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode));
  rmdir(report_directory.c_str());
}

TEST_F(BartTestHelperTest, ReInitalizeReport) {
  test_helper.ReInit(true, gold_files_directory);
  std::string report_directory = test_helper.GetReportDirectory();
  // Verify report directory name
  EXPECT_THAT(report_directory,
              MatchesRegex(gold_files_directory + "........_...._.._report"));
}

TEST_F(BartTestHelperTest, InitalizeBadDirectory) {
  std::string bad_directory = "testing_data/";
  ASSERT_THROW(btest::BartTestHelper new_test_helper(true, bad_directory),
               std::runtime_error);
}

TEST_F(BartTestHelperTest, LogTest) {
  std::string filename = "test_log_file";
  test_helper.OpenLog("test_log_file");
  dealii::deallog << "test" << std::endl;
  EXPECT_TRUE(test_helper.IsLogging());
  test_helper.CloseLog();
  EXPECT_FALSE(test_helper.IsLogging());
  std::ifstream input_log_file(filename);
  std::string line;
  ASSERT_TRUE(input_log_file);
  std::getline(input_log_file, line);
  EXPECT_EQ(line, "DEAL::test");
  remove(filename.c_str());
}

TEST_F(BartTestHelperTest, LogTestThrow) {
  std::string filename = "test_log_file";
  test_helper.OpenLog(filename);
  EXPECT_THROW(test_helper.OpenLog("other_log"), std::runtime_error);
  test_helper.CloseLog();
  remove(filename.c_str());
}

