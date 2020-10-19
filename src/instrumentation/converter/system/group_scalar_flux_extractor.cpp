#include "instrumentation/converter/system/group_scalar_flux_extractor.hpp"

namespace bart::instrumentation::converter::system {

moments::MomentVector GroupScalarFluxExtractor::Convert(
    const moments::SphericalHarmonicI &input) const {
  return input.GetMoment({group_to_extract_, 0, 0});
}

} // namespace bart::instrumentation::converter::system
