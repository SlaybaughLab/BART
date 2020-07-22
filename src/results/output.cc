#include "results/output.h"

#include <deal.II/base/exceptions.h>

namespace bart {

namespace results {

void Output::WriteVector(std::ostream &output_stream,
                         const std::vector<double> to_write) const {
  for (std::vector<double>::size_type i =0 ; i < to_write.size(); ++i)
    output_stream << i << ", " << to_write.at(i) << "\n";
}

void Output::WriteVector(std::ostream &output_stream,
                         std::vector<double> to_write,
                         std::vector<std::string> headers) const {
  AssertThrow(static_cast<int>(headers.size()) == 2,
              dealii::ExcMessage("Error in WriteVector, header size must be "
                                 "size 2"))
  output_stream << headers.at(0) << "," << headers.at(1) << "\n";
  WriteVector(output_stream, to_write);
}

} // namespace results

} // namespace bart
