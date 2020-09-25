#ifndef BART_SRC_INSTRUMENTATION_FACTORY_CONVERTER_FACTORY_REGISTRATION_H_
#define BART_SRC_INSTRUMENTATION_FACTORY_CONVERTER_FACTORY_REGISTRATION_H_

#include "instrumentation/factory/converter_types.h"
#include "instrumentation/factory/converter_factory.h"

namespace bart {

namespace instrumentation {

namespace factory {

namespace registrations {

template <typename T>
class ConverterFactoryRegistration {
 public:
  using InputType = typename T::InputType;
  using OutputType = typename T::OutputType;
  ConverterFactoryRegistration(
      const ConverterName type,
      ConverterInstanceGenerator<typename T::InputType, typename T::OutputType> generator) {
    ConverterFactory<InputType, OutputType>::get().RegisterGenerator(
        type, generator);
  }
};

} // namespace registrations

} // namespace factory

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_FACTORY_CONVERTER_FACTORY_REGISTRATION_H_
