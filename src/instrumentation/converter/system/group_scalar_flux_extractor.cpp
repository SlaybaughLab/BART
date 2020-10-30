#include "instrumentation/converter/system/group_scalar_flux_extractor.hpp"
#include "instrumentation/converter/factory.hpp"

namespace bart::instrumentation::converter::system {

namespace  {
using Factory = ConverterIFactory<moments::SphericalHarmonicI, moments::MomentVector, const int>;
using ConverterType = ConverterI<moments::SphericalHarmonicI, moments::MomentVector>;
auto constructor_function = [](const int group) -> std::unique_ptr<ConverterType> {
  return std::make_unique<GroupScalarFluxExtractor>(group);
};
} // namespace

moments::MomentVector GroupScalarFluxExtractor::Convert(
    const moments::SphericalHarmonicI &input) const {
  return input.GetMoment({group_to_extract_, 0, 0});
}

bool GroupScalarFluxExtractor::is_registered_ = Factory::get()
    .RegisterConstructor(ConverterName::kGroupScalarFluxExtractor,
                         constructor_function);

} // namespace bart::instrumentation::converter::system
