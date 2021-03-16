#ifndef BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_

#include "utility/has_description.h"

namespace bart::calculator::residual {

class DomainIsotropicResidualI : public utility::HasDescription {
 public:
  virtual ~DomainIsotropicResidualI() = default;
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_
