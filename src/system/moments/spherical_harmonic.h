#ifndef BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_H_
#define BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_H_

#include "system/moments/spherical_harmonic_types.h"
#include "system/moments/spherical_harmonic_i.h"

namespace bart {

namespace system {

namespace moments {

class SphericalHarmonic : public SphericalHarmonicI {
 public:
  SphericalHarmonic(const int total_groups,
                    const int max_harmonic_l)
      : total_groups_(total_groups),
        max_harmonic_l_(max_harmonic_l) {};
  virtual ~SphericalHarmonic() = default;

  int total_groups() const override { return total_groups_; }
  int max_harmonic_l() const override { return max_harmonic_l_;}
 private:
  const int total_groups_ = 0;
  const int max_harmonic_l_ = 0;
};

} // namespace moments

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_H_