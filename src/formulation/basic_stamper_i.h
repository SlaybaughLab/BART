#ifndef BART_SRC_BART_FORMULATIOIN_BASIC_STAMPER_H_I_
#define BART_SRC_BART_FORMULATIOIN_BASIC_STAMPER_H_I_

#include "formulation/stamper_i.h"

namespace bart {

namespace formulation {

/*! Basic stamper class for CFEM formulations.
 *
 * The basic stamper stamps the following terms:
 *
 * Bilinear (take a matrix as an input)
 *  - Streaming
 *  - Collision
 *  - Boundary
 *
 * Linear (take a vector as an input)
 *  - Fixed source
 *  - Fission source
 *  - Scattering source
 *
 * Source terms are calculated using moments.
 *
 */
class BasicStamperI : public StamperI {
 public:
  virtual ~BasicStamperI() = default;
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_BART_FORMULATIOIN_BASIC_STAMPER_H_I_