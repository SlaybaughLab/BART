#include "test_helpers/bart_test_helper.h"

#include <exception>
#include <fstream>
#include <sys/stat.h>

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class BARTTestHelperTest : public ::testing::Test {
 protected:
  virtual void TearDown() override;
  std::string gold_files_directory = "test_data/";
  test_helpers::BARTTestHelper test_helper;
};

void BARTTestHelperTest::TearDown() {
    std::string report_directory = test_helper.GetReportDirectory();
    if (report_directory != "")
      rmdir(test_helper.GetReportDirectory().c_str());
}

TEST_F(BARTTestHelperTest, InitializationNoReport) {
  std::string directory_name = test_helper.GetReportDirectory();
  ASSERT_TRUE(directory_name.empty());
}

TEST_F(BARTTestHelperTest, ConstructorReport) {
  test_helpers::BARTTestHelper test_helper2(true, gold_files_directory);
  std::string report_directory = test_helper2.GetReportDirectory();
  // Verify report directory name
  std::string regex = gold_files_directory + "........_...._.._report";
  EXPECT_THAT(report_directory, ::testing::MatchesRegex(regex));
  // Verify report directory existence
  struct stat sb;
  bool status = (stat(report_directory.c_str(), &sb) == 0) && S_ISDIR(sb.st_mode);  
  ASSERT_TRUE(status);
  rmdir(report_directory.c_str());
}

TEST_F(BARTTestHelperTest, ReInitalizeReport) {
  test_helper.ReInit(true, gold_files_directory);
  std::string report_directory = test_helper.GetReportDirectory();
  std::string regex = gold_files_directory + "........_...._.._report";
  // Verify report directory name
  EXPECT_THAT(report_directory, ::testing::MatchesRegex(regex));
}

TEST_F(BARTTestHelperTest, SetReport) {
  ASSERT_EQ(test_helper.GetReportDirectory(), "");
  test_helper.SetReport(true);
  std::string report_directory = test_helper.GetReportDirectory();
  std::string regex = gold_files_directory + "........_...._.._report";
  // Verify report directory name
  EXPECT_THAT(report_directory, ::testing::MatchesRegex(regex));
}

TEST_F(BARTTestHelperTest, SetGoldFilesDirectory) {
  test_helper.ReInit(false, "other_directory/");
  test_helper.SetGoldFilesDirectory(gold_files_directory);
  test_helper.SetReport(true);
  std::string report_directory = test_helper.GetReportDirectory();
  std::string regex = gold_files_directory + "........_...._.._report";
  // Verify report directory name
  EXPECT_THAT(report_directory,::testing::MatchesRegex(regex));
}


TEST_F(BARTTestHelperTest, InitalizeBadDirectory) {
  std::string bad_directory = "testing_data/";
  ASSERT_THROW(test_helpers::BARTTestHelper new_test_helper(true, bad_directory),
               std::runtime_error);
}

TEST_F(BARTTestHelperTest, LogTest) {
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

TEST_F(BARTTestHelperTest, LogTestThrow) {
  std::string filename = "test_log_file";
  test_helper.OpenLog(filename);
  EXPECT_THROW(test_helper.OpenLog("other_log"), std::runtime_error);
  test_helper.CloseLog();
  remove(filename.c_str());
}

} // namespace