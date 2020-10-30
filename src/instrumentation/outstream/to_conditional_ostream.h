#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_CONDITIONAL_OSTREAM_H_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_CONDITIONAL_OSTREAM_H_

#include <deal.II/base/conditional_ostream.h>

#include "instrumentation/outstream/outstream_i.h"
#include "instrumentation/outstream/factory.hpp"

#include "utility/named_type.h"

namespace bart {

namespace instrumentation {

namespace outstream {

class ToConditionalOstream : public OutstreamI<std::string> {
 public:
  using ConditionalOstreamType = dealii::ConditionalOStream;
  using ConditionalOstreamPtrType = std::unique_ptr<ConditionalOstreamType>;
  ToConditionalOstream(ConditionalOstreamPtrType conditional_ostream_ptr)
      : conditional_ostream_ptr_(std::move(conditional_ostream_ptr)) {}
  ToConditionalOstream& Output(const std::string& to_output) override {
    *conditional_ostream_ptr_ << to_output;
    return *this;
  };

  ConditionalOstreamType* conditional_ostream_ptr() const {
    return conditional_ostream_ptr_.get(); }
 private:
  ConditionalOstreamPtrType conditional_ostream_ptr_ = nullptr;
  static bool is_registered_;
};

} // namespace outstream

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_CONDITIONAL_OSTREAM_H_
