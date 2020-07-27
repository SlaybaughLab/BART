#ifndef BART_SRC_INSTRUMENTATION_OUTPUT_TESTS_OUTPUT_MOCK_H_
#define BART_SRC_INSTRUMENTATION_OUTPUT_TESTS_OUTPUT_MOCK_H_

#include "instrumentation/output/output_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace instrumentation {

namespace output {

template <typename DataType>
class OutputMock : public OutputI<DataType> {
 public:
  MOCK_METHOD(OutputMock&, Output, (const DataType&), (override));
};

} // namespace output

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTPUT_TESTS_OUTPUT_MOCK_H_
