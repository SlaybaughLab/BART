#ifndef BART_SRC_INSTRUMENTATION_BASIC_INSTRUMENT_H_
#define BART_SRC_INSTRUMENTATION_BASIC_INSTRUMENT_H_

#include <memory>

#include "instrument_i.h"
#include "instrumentation/output/output_i.h"
#include "utility/has_dependencies.h"
#include "utility/has_description.h"

namespace bart {

namespace instrumentation {

/* \brief A basic instrument has the same data type for both reading
 * and outputting and no converter. */
template <typename InputType>
class BasicInstrument : public InstrumentI<InputType>,
                        public utility::HasDependencies,
                        public utility::HasDescription {
 public:
  using OutputterType = typename output::OutputI<InputType>;
  explicit BasicInstrument(std::unique_ptr<OutputterType> outputter_ptr)
      : outputter_ptr_(std::move(outputter_ptr)) {
    AssertPointerNotNull(outputter_ptr_.get(), __func__ , "outputter");
    set_description("Basic instrument", utility::DefaultImplementation(true));
  }
  void Read(const InputType &input) override {}
  virtual OutputterType* outputter_ptr() { return outputter_ptr_.get(); }
 protected:
  std::unique_ptr<OutputterType> outputter_ptr_ = nullptr;
};

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_BASIC_INSTRUMENT_H_
