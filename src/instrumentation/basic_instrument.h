#ifndef BART_SRC_INSTRUMENTATION_BASIC_INSTRUMENT_H_
#define BART_SRC_INSTRUMENTATION_BASIC_INSTRUMENT_H_

#include <memory>
#include <string>

#include "instrument_i.h"
#include "instrumentation/outstream/outstream_i.h"
#include "utility/has_dependencies.h"
#include "utility/has_description.h"

namespace bart {

namespace instrumentation {

/* \brief A basic instrument has the same data type for both reading
 * and outputting and no converter. */
template <typename InputType = std::string>
class BasicInstrument : public InstrumentI<InputType>,
                        public utility::HasDependencies,
                        public utility::HasDescription {
 public:
  using OutstreamType = typename outstream::OutstreamI<InputType>;
  explicit BasicInstrument(std::unique_ptr<OutstreamType> outstream_ptr)
      : outstream_ptr_(std::move(outstream_ptr)) {
    AssertPointerNotNull(outstream_ptr_.get(), "outstream_ptr", __func__);
    set_description("Basic instrument", utility::DefaultImplementation(true));
  }
  void Read(const InputType &input) override { outstream_ptr_->Output(input); }
  virtual OutstreamType* outstream_ptr() { return outstream_ptr_.get(); }
 protected:
  std::unique_ptr<OutstreamType> outstream_ptr_ = nullptr;
};

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_BASIC_INSTRUMENT_H_
