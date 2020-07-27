#ifndef BART_SRC_INSTRUMENTATION_OUTPUT_OUTPUT_I_H_
#define BART_SRC_INSTRUMENTATION_OUTPUT_OUTPUT_I_H_

namespace bart {

namespace instrumentation {

namespace output {

template <typename DataType>
class OutputI {
 public:
  virtual ~OutputI() = default;
  virtual OutputI& Output(const DataType& to_output) = 0;
};

} // namespace output

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_OUTPUT_OUTPUT_I_H_
