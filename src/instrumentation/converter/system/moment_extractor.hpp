#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_SYSTEM_MOMENT_EXTRACTOR_HPP_
#define BART_SRC_INSTRUMENTATION_CONVERTER_SYSTEM_MOMENT_EXTRACTOR_HPP_

#include "instrumentation/converter/converter_i.h"

#include "system/moments/spherical_harmonic_i.h"

namespace bart::instrumentation::converter::system {

namespace moments = bart::system::moments;

/*! \brief Extracts a moment vector from a collection. */
class MomentExtractor :
    public ConverterI<moments::SphericalHarmonicI,
                      moments::MomentVector> {
 public:
  moments::MomentVector Convert(const moments::SphericalHarmonicI &input) const override;
};

} // namespace bart::instrumentation::converter::system

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_SYSTEM_MOMENT_EXTRACTOR_HPP_
