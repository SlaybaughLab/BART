#include "instrumentation/outstream/vector_to_vtu.hpp"

namespace bart::instrumentation::outstream {

auto VectorToVTU::Output(const Vector& to_output) -> VectorToVTU& {
  return *this;
}
} // namespace bart::instrumentation::outstream

