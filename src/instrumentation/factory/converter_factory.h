#ifndef BART_SRC_INSTRUMENTATION_FACTORY_CONVERTER_FACTORY_H_
#define BART_SRC_INSTRUMENTATION_FACTORY_CONVERTER_FACTORY_H_

#include <unordered_map>
#include <memory>

#include "instrumentation/factory/converter_types.h"

namespace bart {

namespace instrumentation {

namespace converter {
template <typename InputType, typename OutputType>
class ConverterI;
} // namespace converter

namespace factory {

template <typename InputType, typename OutputType>
using ConverterInstanceGenerator =
    std::unique_ptr<converter::ConverterI<InputType, OutputType>>(*)();

template <typename InputType, typename OutputType>
class ConverterFactory {
 public:
  static ConverterFactory& get() {
    static ConverterFactory instance;
    return instance;
  };
  ConverterInstanceGenerator<InputType, OutputType> get_generator(
      const ConverterName type) {
    return converter_generators_.at(type);
  }

  bool RegisterGenerator(
      const ConverterName type,
      const ConverterInstanceGenerator<InputType, OutputType>& generator) {
    return converter_generators_.insert(std::make_pair(type, generator)).second;
  }

 private:
  ConverterFactory() = default;
  ConverterFactory(const ConverterFactory&);
  ~ConverterFactory() = default;

  std::unordered_map<ConverterName,
                     ConverterInstanceGenerator<InputType, OutputType>>
      converter_generators_;

};

} // namespace factory

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_FACTORY_CONVERTER_FACTORY_H_
