#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_OUTSTREAM_I_H_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_OUTSTREAM_I_H_

namespace bart {

namespace instrumentation {

namespace outstream {

template <typename DataType>
class OutstreamI {
 public:
  ~OutstreamI() = default;
  virtual OutstreamI& Output(const DataType& to_output) = 0;
};

} // namespace outstream

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_OUTSTREAM_I_H_
