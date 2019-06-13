#include "system/moments/spherical_harmonic.h"

namespace bart {

namespace system {

namespace moments {

SphericalHarmonic::SphericalHarmonic(const int total_groups,
                                     const int max_harmonic_l)
    : total_groups_(total_groups),
      max_harmonic_l_(max_harmonic_l) {
  for (int group = 0; group < total_groups; ++group) {
    for (int l = 0; l <= max_harmonic_l; ++l) {
      for (int m = -l; m <= l; ++m) {
        moments_.insert(
            std::make_pair<MomentIndex, MomentVector>({group, l, m}, {}));
      }
    }
  }
}
} // namespace moments

} // namespace system

} // namespace bart