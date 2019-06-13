#ifndef BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_
#define BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_

#include "system/moments/spherical_harmonic_types.h"
#include "system/moments/spherical_harmonic_types.h"

namespace bart {

namespace system {

namespace moments {

/*! \brief Interface for spherical harmonic moments storage class
 *
 */
class SphericalHarmonicI {
 public:
  virtual ~SphericalHarmonicI() = default;

  //virtual MomentVector& operator[](const MomentIndex index) = 0;

  virtual int total_groups() const = 0;
  virtual int max_harmonic_l() const = 0;
};

} // namespace moments

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_