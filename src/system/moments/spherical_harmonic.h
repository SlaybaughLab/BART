#ifndef BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_H_
#define BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_H_

#include "system/moments/spherical_harmonic_types.h"
#include "system/moments/spherical_harmonic_i.h"

namespace bart {

namespace system {

namespace moments {

class SphericalHarmonic : public SphericalHarmonicI {
 public:
  virtual ~SphericalHarmonic() = default;
};

} // namespace moments

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_H_