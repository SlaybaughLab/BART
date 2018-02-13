#include "gtest/gtest.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/logstream.h>

#include <iostream>
#include <memory>
#include <fstream>

#include "../log_handler.h"

class LogHandlerTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    log_stream = std::make_unique<std::ofstream>("log_file");
  }
  btest::LogHandler test_handler;
  std::unique_ptr<std::ofstream> log_stream;
};

TEST_F(LogHandlerTest, AttachToDealLog) {
  test_handler.OpenLogFile(std::move(log_stream));
  dealii::deallog << "test" << std::endl;
  EXPECT_TRUE(test_handler.IsLogging());
  test_handler.CloseLogFile();
  EXPECT_FALSE(test_handler.IsLogging());
  std::ifstream input_log_file("log_file");
  std::string line;
  ASSERT_TRUE(input_log_file);
  std::getline(input_log_file, line);
  EXPECT_EQ(line, "DEAL::test");
  remove("log_file");
}
