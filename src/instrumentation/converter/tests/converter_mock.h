#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_TESTS_CONVERTER_MOCK_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_TESTS_CONVERTER_MOCK_H_

#include "instrumentation/converter/converter_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace instrumentation {

namespace converter {

template <typename InputType, typename OutputType>
class ConverterMock : public ConverterI<InputType, OutputType> {
 public:
  MOCK_METHOD(OutputType, Convert, (const InputType&), (const, override));
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_TESTS_CONVERTER_MOCK_H_
