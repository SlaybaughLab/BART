#include "utility/reporter/mpi.h"

#include <deal.II/base/conditional_ostream.h>
#include <gtest/gtest.h>

#include "convergence/status.h"
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

} // namespace