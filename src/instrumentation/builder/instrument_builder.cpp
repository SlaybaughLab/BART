#include "instrumentation/builder/instrument_builder.hpp"

#include <deal.II/base/mpi.h>
#include <deal.II/base/conditional_ostream.h>

#include "instrumentation/basic_instrument.h"
#include "instrumentation/instrument.h"

#include "instrumentation/converter/factory.h"
#include "instrumentation/outstream/factory.h"
#include "utility/colors.h"

namespace bart::convergence {
struct Status;
} // namespace bart::convergence

namespace bart::instrumentation::builder {

namespace  {

using ConditionalOstream = dealii::ConditionalOStream;
using ConditionalOstreamPtr = std::unique_ptr<ConditionalOstream>;
using ConvergenceStatus = convergence::Status;
using ConvergenceInstrument = InstrumentI<convergence::Status>;
using ConverterName = instrumentation::converter::ConverterName;
using StringInstrument = InstrumentI<std::string>;
using StringColorPair = std::pair<std::string, utility::Color>;

auto GetConditionalOstream = []() {
  return outstream::OutstreamIFactory<std::string, ConditionalOstreamPtr>::get()
    .GetConstructor(outstream::OutstreamName::kToConditionalOstream)(
        std::make_unique<ConditionalOstream>(
            std::cout,
            dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));
};

} // namespace

// Convergence Status

template <>
auto InstrumentBuilder::BuildInstrument<ConvergenceStatus>(
    const InstrumentName name) -> std::unique_ptr<ConvergenceInstrument> {
  switch (name) {
    case InstrumentName::kConvergenceStatusToConditionalOstream: {
      return std::make_unique<Instrument<ConvergenceStatus, std::string>>(
          converter::ConverterIFactory<ConvergenceStatus, std::string>::get()
              .GetConstructor(ConverterName::kConvergenceToString)(),
          GetConditionalOstream());
    }
    default:
    AssertThrow(false,
                dealii::ExcMessage("Bad instrument name passed to builder"))
  }
}

// STRING-COLOR-PAIR ===========================================================

template<>
std::unique_ptr<InstrumentI<StringColorPair>> InstrumentBuilder::BuildInstrument(
    const InstrumentName name) {

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
auto InstrumentBuilder::BuildInstrument(const InstrumentName name)
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
