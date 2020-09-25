#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_CONVERTER_I_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_CONVERTER_I_H_

#include "utility/has_description.h"
#include "instrumentation/factory/converter_factory_registration.h"

namespace bart {

namespace instrumentation {

namespace converter {

/* \brief Interface for classes that convert one datatype to another.
 *
 * Classes that implement this interface are designed to work within an
 * instrument to convert an incoming data type to the type required by the
 * outputter
 */
template <typename InputType, typename OutputType>
class ConverterI {
 public:
  virtual ~ConverterI() = default;
  virtual OutputType Convert(const InputType& input) const = 0;
};

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_CONVERTER_I_H_
