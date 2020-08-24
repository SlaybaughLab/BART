#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_TESTS_OUTSTREAM_MOCK_H_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_TESTS_OUTSTREAM_MOCK_H_

#include "instrumentation/outstream/outstream_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace instrumentation {

namespace outstream {

template <typename DataType>
class OutstreamMock : public OutstreamI<DataType> {
 public:
  MOCK_METHOD(OutstreamMock&, Output, (const DataType&), (override));
};

} // namespace outstream

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_TESTS_OUTSTREAM_MOCK_H_
