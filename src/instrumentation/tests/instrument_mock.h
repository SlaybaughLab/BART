#ifndef BART_SRC_INSTRUMENTATION_TESTS_INSTRUMENT_MOCK_H_
#define BART_SRC_INSTRUMENTATION_TESTS_INSTRUMENT_MOCK_H_

#include "instrumentation/instrument_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace instrumentation {

template <typename InputType>
class InstrumentMock : public InstrumentI<InputType> {
 public:
  MOCK_METHOD(void, Read, (const InputType& input), (override));
};

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_TESTS_INSTRUMENT_MOCK_H_
