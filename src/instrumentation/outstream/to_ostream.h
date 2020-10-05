#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_OSTREAM_H_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_OSTREAM_H_

#include "instrumentation/outstream/outstream_i.h"
#include "instrumentation/outstream/factory.h"

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
  inline static bool is_registered_ =
      OutstreamIFactory<std::string, std::unique_ptr<std::ostream>>::get()
      .RegisterConstructor(OutstreamName::kToOstream,
          [](std::unique_ptr<std::ostream> ostream_ptr) {
        std::unique_ptr<OutstreamI<std::string>> return_ptr =
            std::make_unique<ToOstream>(std::move(ostream_ptr));
        return return_ptr;
          });
};

} // namespace outstream

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_TO_OSTREAM_H_
