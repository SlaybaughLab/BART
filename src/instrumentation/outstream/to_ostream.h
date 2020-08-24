#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_OSTREAM_H_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_OSTREAM_H_

#include "instrumentation/outstream/outstream_i.h"

#include <memory>
#include <ostream>

namespace bart {

namespace instrumentation {

namespace outstream {

class ToOstream : public OutstreamI<std::string> {
 public:
  ToOstream(std::unique_ptr<std::ostream> ostream_ptr)
      : ostream_ptr_(std::move(ostream_ptr)) {}
  ToOstream& Output(const std::string &to_output) override {
    *ostream_ptr_ << to_output;
    return *this;
  }

  std::ostream* ostream_ptr() { return ostream_ptr_.get(); }

 private:
  std::unique_ptr<std::ostream> ostream_ptr_;
};

} // namespace outstream

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_OSTREAM_H_
