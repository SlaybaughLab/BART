#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_TO_VTU_HPP_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_TO_VTU_HPP_

#include <deal.II/lac/vector.h>

#include "instrumentation/outstream/outstream_i.h"

namespace bart::instrumentation::outstream {

class VectorToVTU : public OutstreamI<dealii::Vector<double>> {
 public:
  using Vector = dealii::Vector<double>;
  auto Output(const Vector& to_output) -> VectorToVTU& override;
};

} // namespace bart::instrumentation::outstream

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_VECTOR_TO_VTU_HPP_
