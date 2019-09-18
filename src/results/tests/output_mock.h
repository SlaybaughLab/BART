#ifndef BART_SRC_RESULTS_TESTS_OUTPUT_MOCK_H_
#define BART_SRC_RESULTS_TESTS_OUTPUT_MOCK_H_

#include "results/output_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace results {


class OutputMock : public OutputI {
 public:
  MOCK_METHOD1(AddData, void(system::System&));
  MOCK_CONST_METHOD1(WriteData, void(std::ostream&));
};

} // namespace results

} // namespace bart

#endif //BART_SRC_RESULTS_TESTS_OUTPUT_MOCK_H_
