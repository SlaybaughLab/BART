#include "utility/reporter/mpi.h"

#include <deal.II/base/conditional_ostream.h>
#include <gtest/gtest.h>

#include "utility/reporter/colors.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class UtilityReporterMpiTest : public ::testing::Test {
 protected:
  std::ostringstream string_stream;
  std::unique_ptr<utility::reporter::Mpi> test_reporter_;
  void SetUp() override;
};

void UtilityReporterMpiTest::SetUp() {
  auto pout_ptr = std::make_unique<dealii::ConditionalOStream>(string_stream, true);
  test_reporter_ = std::make_unique<utility::reporter::Mpi>(std::move(pout_ptr));
}

/* Verifies that the reporter sends the string to the ostream */
TEST_F(UtilityReporterMpiTest, StringReport) {
  std::string to_report = "report me";
  test_reporter_->Report(to_report);
  EXPECT_EQ(string_stream.str(), to_report);
}

TEST_F(UtilityReporterMpiTest, InstreamOperator) {
  std::string to_report = "report me";
  *test_reporter_ << to_report;
  EXPECT_EQ(string_stream.str(), to_report);
}

TEST_F(UtilityReporterMpiTest, ColoredStringReport) {
  using Color = utility::reporter::Color;
  std::unordered_map<Color, std::string> color_string{
      {Color::Reset, "\033[0m"},
      {Color::Red,   "\033[31m"},
      {Color::Green, "\033[32m"},
      {Color::Blue,  "\033[34m"},
  };

  for (auto color_pair : color_string) {
    string_stream.str("");
    std::string to_report = "report me";
    std::string colored_report = color_pair.second + to_report + color_string.at(Color::Reset);
    test_reporter_->Report(to_report, color_pair.first);
    EXPECT_EQ(string_stream.str(), colored_report);
  }

}

TEST_F(UtilityReporterMpiTest, InstreamOperatorColor) {
  using Color = utility::reporter::Color;
  std::unordered_map<Color, std::string> color_string{
      {Color::Reset, "\033[0m"},
      {Color::Red,   "\033[31m"},
      {Color::Green, "\033[32m"},
      {Color::Blue,  "\033[34m"},
  };

  for (auto color_pair : color_string) {
    string_stream.str("");
    std::string to_report = "report me";
    std::string colored_report = color_pair.second + to_report + color_string.at(Color::Reset);
    *test_reporter_ << color_pair.second << to_report << Color::Reset;
    EXPECT_EQ(string_stream.str(), colored_report);
  }

}

} // namespace