#include "instrumentation/converter/system/moment_extractor.hpp"

namespace bart::instrumentation::converter::system {

moments::MomentVector MomentExtractor::Convert(const moments::SphericalHarmonicI &input) const {
  return bart::system::moments::MomentVector();
}
} // namespace bart::instrumentation::converter::system
