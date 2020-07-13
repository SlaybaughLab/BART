#ifndef BART_SRC_RESULTS_TESTS_OUTPUT_MOCK_H_
#define BART_SRC_RESULTS_TESTS_OUTPUT_MOCK_H_

#include "results/output_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace results {


class OutputMock : public OutputI {
 public:
  MOCK_METHOD(void, AddData, (system::System&), (override));
  MOCK_METHOD(void, WriteData, (std::ostream&), (override, const));
  MOCK_METHOD(void, WriteMasterFile, (std::ostream &output_stream,
      std::vector<std::string> filenames), (override, const));
  MOCK_METHOD(void, WriteVector, (std::ostream &, const std::vector<double>),
              (const, override));
};

} // namespace results

} // namespace bart

#endif //BART_SRC_RESULTS_TESTS_OUTPUT_MOCK_H_
