#include "prm_parameters.h"

namespace bart {

namespace problem {

PrmParameters::PrmParameters(const dealii::ParameterHandler &handler) {
  dimension_ = handler.get_integer("problem dimension");
}

} // namespace problem

} // namespace bart
