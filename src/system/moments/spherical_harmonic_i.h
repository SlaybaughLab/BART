#ifndef BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_
#define BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_

#include "system/moments/moment_types.h"

namespace bart {

namespace system {

namespace moments {

/*! \brief Interface for spherical harmonic moments storage class
 *
 */
class SphericalHarmonicI {
 public:
  virtual ~SphericalHarmonicI() = default;
};

} // namespace moments

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_