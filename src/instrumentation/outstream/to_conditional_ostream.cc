#include "instrumentation/outstream/to_conditional_ostream.h"

namespace bart::instrumentation::outstream {

bool ToConditionalOstream::is_registered_ =
    OutstreamIFactory<std::string, ConditionalOstreamPtrType>::get()
        .RegisterConstructor(OutstreamName::kToConditionalOstream,
                             [](ConditionalOstreamPtrType conditional_ostream_ptr)
                                 -> std::unique_ptr<OutstreamI<std::string>> {
                               return std::make_unique<ToConditionalOstream>(
                                   std::move(conditional_ostream_ptr));
                             });

} // namespace bart::instrumentation::outstream
