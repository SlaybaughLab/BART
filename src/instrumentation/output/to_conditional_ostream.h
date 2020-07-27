#ifndef BART_SRC_INSTRUMENTATION_OUTPUT_TO_CONDITIONAL_OSTREAM_H_
#define BART_SRC_INSTRUMENTATION_OUTPUT_TO_CONDITIONAL_OSTREAM_H_

#include "instrumentation/output/output_i.h"

#include <deal.II/base/conditional_ostream.h>

namespace bart {

namespace instrumentation {

namespace output {

class ToConditionalOstream : public OutputI<std::string> {
 public:
  using ConditionalOstreamType = dealii::ConditionalOStream;
  using ConditionalOstreamPtrType = std::unique_ptr<ConditionalOstreamType>;
  ToConditionalOstream(ConditionalOstreamPtrType conditional_ostream_ptr)
      : conditional_ostream_ptr_(std::move(conditional_ostream_ptr)) {}
  void Output(const std::string& to_output) override {};

  ConditionalOstreamType* conditional_ostream_ptr() const {
    return conditional_ostream_ptr_.get(); }
 private:
  ConditionalOstreamPtrType conditional_ostream_ptr_ = nullptr;

};

} // namespace output

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTPUT_TO_CONDITIONAL_OSTREAM_H_
