#ifndef BART_SRC_INSTRUMENTATION_INSTRUMENT_I_H_
#define BART_SRC_INSTRUMENTATION_INSTRUMENT_I_H_

namespace bart {

namespace instrumentation {

template <typename InputType>
class InstrumentI {
 public:
  virtual ~InstrumentI() = default;
  virtual void Read(const InputType& input) = 0;
};

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_INSTRUMENT_I_H_
