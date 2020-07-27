#ifndef BART_SRC_INSTRUMENTATION_OUTPUT_TO_CONDITIONAL_OSTREAM_H_
#define BART_SRC_INSTRUMENTATION_OUTPUT_TO_CONDITIONAL_OSTREAM_H_

#include "instrumentation/output/output_i.h"

#include <deal.II/base/conditional_ostream.h>

namespace bart {

namespace instrumentation {

namespace output {

template <typename DataType>
class ToConditionalOstream : public OutputI<DataType> {
 public:
  void Output(const DataType& to_output) override {};
};

} // namespace output

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTPUT_TO_CONDITIONAL_OSTREAM_H_
