#include "instrumentation/outstream/to_ostream.hpp"

namespace bart::instrumentation::outstream {

bool ToOstream::is_registered_ =
    OutstreamIFactory<std::string, std::unique_ptr<std::ostream>>::get()
        .RegisterConstructor(OutstreamName::kToOstream,
                             [](std::unique_ptr<std::ostream> ostream_ptr) {
                               std::unique_ptr<OutstreamI<std::string>> return_ptr =
                                   std::make_unique<ToOstream>(std::move(ostream_ptr));
                               return return_ptr;
                             });

} // namespace bart::instrumentation::outstream
