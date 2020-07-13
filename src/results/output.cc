#include "results/output.h"

namespace bart {

namespace results {

void Output::WriteVector(const std::vector<double> to_write,
    std::ostream &output_stream) const {
  for (std::vector<double>::size_type i =0 ; i < to_write.size(); ++i)
    output_stream << i << ", " << to_write.at(i) << "\n";
}

} // namespace results

} // namespace bart
