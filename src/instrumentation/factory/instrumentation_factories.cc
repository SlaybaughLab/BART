#include "instrumentation_factories.h"

#include <deal.II/base/conditional_ostream.h>

#include "instrumentation/output/to_conditional_ostream.h"

namespace bart {

namespace instrumentation {

namespace factory {

template <>
std::unique_ptr<instrumentation::output::OutputI<std::string>>
MakeOutputter<std::string, std::unique_ptr<dealii::ConditionalOStream>>(
    std::unique_ptr<dealii::ConditionalOStream> conditional_ostream_ptr) {
  return std::make_unique<instrumentation::output::ToConditionalOstream>(
      std::move(conditional_ostream_ptr));
}

} // namespace factory

} // namespace instrumentation

} // namespace bart
