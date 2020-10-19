#include "instrumentation/converter/system/moment_extractor.hpp"

namespace bart::instrumentation::converter::system {

moments::MomentVector MomentExtractor::Convert(
    const moments::SphericalHarmonicI &input) const {
  return input.GetMoment({group_to_extract_, 0, 0});
}

} // namespace bart::instrumentation::converter::system
