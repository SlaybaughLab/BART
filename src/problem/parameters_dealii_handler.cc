#include "parameters_dealii_handler.h"

namespace bart {

namespace problem {

ParametersDealiiHandler::ParametersDealiiHandler(
    const dealii::ParameterHandler &handler) {
  spatial_dimension_ = handler.get_integer("problem dimension");
}

} // namespace problem

} // namespace bart
