#include "instrumentation/instrument_array.hpp"

#include "system/moments/spherical_harmonic_i.h"

namespace bart::instrumentation {

template<typename InputType>
void InstrumentArray<InputType>::Read(const InputType &/*input*/) {}


template<typename InputType>
InstrumentArray<InputType>& InstrumentArray<InputType>::AddInstrument(
    std::unique_ptr<InstrumentType> to_add) {
  instruments_.emplace_back(std::move(to_add));
  return *this;
}

template class InstrumentArray<system::moments::SphericalHarmonicI>;

} // namespace bart::instrumentation
