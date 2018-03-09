#include <gmock/gmock.h>
#include "../bart_test_helper.h"

#include <sys/stat.h>
#include <exception>
#include <fstream>

#include <gtest/gtest.h>

class TestHelperIntTest : public ::testing::Test {
 protected:
  std::string gold_files_directory = "test_data/";
  btest::BARTTestHelper test_helper;
  void TearDown() override;
};

void TestHelperIntTest::TearDown() {
  std::string report_directory = test_helper.GetReportDirectory();
  if (report_directory.empty() &&
      (btest::GlobalBARTTestHelper().GetReportDirectory() != report_directory))
    rmdir(test_helper.GetReportDirectory().c_str());
}

TEST_F(TestHelperIntTest, IntegrationTestGoodNoReport) {
  //btest::BARTTestHelper test_helper(false, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper";
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "1234\n";
  actual_stream.close();
  EXPECT_TRUE(test_helper.GoldTest(filename)) << test_helper.GetFailMessage();
  // Check deleted
  struct stat sb;
  EXPECT_FALSE(stat(filename.c_str(), &sb) == 0);
}

TEST_F(TestHelperIntTest, IntegrationTestBadNoReport) {
  //btest::BARTTestHelper test_helper(false, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper";
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "1235";
  actual_stream.close();
  EXPECT_FALSE(test_helper.GoldTest(filename)) << test_helper.GetFailMessage();
  // Check deleted
  struct stat sb;
  EXPECT_FALSE(stat(filename.c_str(), &sb) == 0);
}

TEST_F(TestHelperIntTest, IntegrationTestNoActual) {
  //btest::BARTTestHelper test_helper(false, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper";
  EXPECT_FALSE(test_helper.GoldTest(filename));
  EXPECT_EQ(test_helper.GetFailMessage(), "Bad test file");
}

TEST_F(TestHelperIntTest, IntegrationTestBadGoldFileNoReport) {
  std::string bad_gold_dir = "/testing_files";
  test_helper.ReInit(false, bad_gold_dir);
  // Make actual file
  std::string filename = "bart_test_helper";
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "1235";
  actual_stream.close();
  EXPECT_FALSE(test_helper.GoldTest(filename));
  EXPECT_EQ(test_helper.GetFailMessage(), "Bad gold file");
  // Check deleted
  struct stat sb;
  EXPECT_FALSE(stat(filename.c_str(), &sb) == 0);
}

TEST_F(TestHelperIntTest, IntegrationTestBadGoldReport) {
  test_helper.ReInit(true, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper_2";
  std::string report_directory = test_helper.GetReportDirectory();
  
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "1235";
  actual_stream.close();
  EXPECT_FALSE(test_helper.GoldTest(filename)) << test_helper.GetFailMessage();
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
  //rmdir(report_directory.c_str());
}

TEST_F(TestHelperIntTest, IntegrationTestBadReport) {
  test_helper.ReInit(true, gold_files_directory);
  // Make actual file
  std::string filename = "bart_test_helper";
  std::string report_directory = test_helper.GetReportDirectory();
  
  std::ofstream actual_stream(filename, std::ios_base::out);
  actual_stream << "12\n34";
  actual_stream.close();
  EXPECT_FALSE(test_helper.GoldTest(filename))  << test_helper.GetFailMessage();
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
  //rmdir(report_directory.c_str());
}  
