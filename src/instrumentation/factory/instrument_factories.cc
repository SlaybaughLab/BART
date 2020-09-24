#include "instrumentation/factory/instrument_factories.h"

namespace bart {

namespace instrumentation {

namespace factory {

std::unique_ptr<FourierInstrumentType> MakeErrorFourierTransformInstrument(
    DealiiVector solution) {
  return std::unique_ptr<FourierInstrumentType>();
}

} // namespace factory

} // namespace instrumentation

} // namespace bart

