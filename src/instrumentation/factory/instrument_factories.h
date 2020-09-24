#ifndef BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENT_FACTORIES_H_
#define BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENT_FACTORIES_H_

#include <deal.II/lac/vector.h>

#include "instrumentation/instrument_i.h"

namespace bart {

namespace instrumentation {

namespace factory {

//using AbsoluteValue =
using DealiiVector = dealii::Vector<double>;
using FourierInstrumentType = instrumentation::InstrumentI<DealiiVector>;

std::unique_ptr<FourierInstrumentType> MakeErrorFourierTransformInstrument(
    DealiiVector solution);

} // namespace factory

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENT_FACTORIES_H_
