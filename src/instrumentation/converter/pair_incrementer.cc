#include "instrumentation/converter/pair_incrementer.h"

#include <complex>
#include <vector>

#include "instrumentation/converter/factory.hpp"

namespace bart {

namespace instrumentation {

namespace converter {

template<typename InputType>
std::pair<int, InputType> PairIncrementer<InputType>::Convert(const InputType &input) const {
  return std::pair<int, InputType>{++increment_, input};
}

template <typename InputType>
bool PairIncrementer<InputType>::is_registered_ =
    converter::ConverterIFactory<InputType, std::pair<int, InputType>>::get()
        .RegisterConstructor(converter::ConverterName::kPairIncrementer,
                             [](){
                               std::unique_ptr<ConverterI<InputType, std::pair<int, InputType>>>
                                   return_ptr = std::make_unique<PairIncrementer<InputType>>();
                               return return_ptr;
                             });

template class PairIncrementer<double>;
template class PairIncrementer<std::string>;
template class PairIncrementer<std::vector<std::complex<double>>>;

} // namespace converter

} // namespace instrumentation

} // namespace bart
