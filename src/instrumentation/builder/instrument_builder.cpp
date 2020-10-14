#include "instrumentation/builder/instrument_builder.hpp"

#include <deal.II/base/mpi.h>
#include <deal.II/base/conditional_ostream.h>

#include "instrumentation/basic_instrument.h"
#include "instrumentation/instrument.h"

#include "instrumentation/converter/factory.h"
#include "instrumentation/outstream/factory.h"
#include "utility/colors.h"

namespace bart::instrumentation::builder {

namespace  {

using ConditionalOstream = dealii::ConditionalOStream;
using ConditionalOstreamPtr = std::unique_ptr<ConditionalOstream>;
using StringInstrument = InstrumentI<std::string>;
using StringColorPair = std::pair<std::string, utility::Color>;

} // namespace

// STRING-COLOR-PAIR ===========================================================

template<>
std::unique_ptr<InstrumentI<StringColorPair>> InstrumentBuilder::BuildInstrument(
    InstrumentName name) {

  switch (name) {
    case InstrumentName::kColorStatusToConditionalOstream: {
      return std::make_unique<Instrument<StringColorPair, std::string>>(
          converter::ConverterIFactory<StringColorPair, std::string>::get()
              .GetConstructor(converter::ConverterName::kStringColorPairToString)(),
          outstream::OutstreamIFactory<std::string, ConditionalOstreamPtr>::get()
              .GetConstructor(outstream::OutstreamName::kToConditionalOstream)(
                  std::make_unique<ConditionalOstream>(
                      std::cout,
                      dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)));
    }
    default:
      AssertThrow(false,
                  dealii::ExcMessage("Bad instrument name passed to builder"))
  }
}

// STRING ======================================================================
template <>
auto InstrumentBuilder::BuildInstrument(InstrumentName name)
-> std::unique_ptr<StringInstrument> {
  switch (name) {
    case InstrumentName::kStringToConditionalOstream: {
      return std::make_unique<BasicInstrument<std::string>>(
          outstream::OutstreamIFactory<std::string, ConditionalOstreamPtr>::get()
              .GetConstructor(outstream::OutstreamName::kToConditionalOstream)(
                  std::make_unique<ConditionalOstream>(
                      std::cout,
                      dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)));
    };
    default:
    AssertThrow(false,
                dealii::ExcMessage("Bad instrument name passed to builder"))
  }
}

} // namespace bart::instrumentation::builder
