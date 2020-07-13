#ifndef BART_SRC_RESULTS_TESTS_OUTPUT_TEST_H_
#define BART_SRC_RESULTS_TESTS_OUTPUT_TEST_H_

#include <ostream>

#include "results/output.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace bart {

namespace results {

namespace testing {

/* This class includes tests for the abstract results::Output class. These can
 * be run by derived classes to verify that the underlying base class
 * functionality is working properly. */
class OutputTest : public ::testing::Test {
 public:
  virtual ~OutputTest() = default;
  // Test writing of vector to a .csv file
  void TestWriteVector(results::Output* test_output);
  void TestWriteVectorWithHeaders(results::Output* test_output);
};

/* Test for TestWriteVector. Verifies that a vector is converted into a csv
 * format properly. We will use a regex to ignore the whitespace and just
 * verify we have each integer index with a comma, followed by the value and
 * a newline. */
void OutputTest::TestWriteVector(results::Output *test_output) {
  std::ostringstream output_string_stream, regex;
  const auto output_vector = test_helpers::RandomVector(4, -100, 100);
  // Create expected string
  for (std::vector<double>::size_type i = 0; i < output_vector.size(); ++i) {
    regex << "\\s*" << i << "\\s*\\,\\s*" << output_vector.at(i) << "\\s*\n\\s*";
  }
  test_output->WriteVector(output_string_stream, output_vector);
  EXPECT_THAT(output_string_stream.str(), ::testing::MatchesRegex(regex.str()))
            << "TestWriteVector failed";
}

/* Test for TestWriteVector with headers.  */
void OutputTest::TestWriteVectorWithHeaders(results::Output *test_output) {
  std::ostringstream output_string_stream, regex;
  const auto output_vector = test_helpers::RandomVector(4, -100, 100);
  const std::vector<std::string> headers{"index", "values"},
      short_header{"index"}, long_header{"index", "values", "extra"};
  // Create expected string
  regex << "\\s*index\\s*\\,\\s*values\\s*\n\\s*";
  for (std::vector<double>::size_type i = 0; i < output_vector.size(); ++i) {
    regex << "\\s*" << i << "\\s*\\,\\s*" << output_vector.at(i) << "\\s*\n\\s*";
  }
  // verify bad header length throws
  for (auto& bad_header : {short_header, long_header}) {
    EXPECT_ANY_THROW(test_output->WriteVector(output_string_stream, output_vector, bad_header))
              << "TestWriteVectorWithHeaders failed bad header throw";
    output_string_stream.str() = "";
  }

  test_output->WriteVector(output_string_stream, output_vector, headers);
  EXPECT_THAT(output_string_stream.str(), ::testing::MatchesRegex(regex.str()))
            << "TestWriteVectorWithHeaders failed regex match";
}

} // namespace testing

} // namespace results

} // namespace bart

#endif //BART_SRC_RESULTS_TESTS_OUTPUT_TEST_H_
