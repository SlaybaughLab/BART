#include "system/moments/spherical_harmonic.hpp"

namespace bart {

namespace system {

namespace moments {

SphericalHarmonic::SphericalHarmonic(const int total_groups,
                                     const int max_harmonic_l)
    : total_groups_(total_groups),
      max_harmonic_l_(max_harmonic_l) {

  AssertThrow(total_groups > 0,
              dealii::ExcMessage("Error system::moments::SphericalHarmonic "
                                 "constructor, total_groups must be > 0"));
  AssertThrow(max_harmonic_l >= 0 ,
              dealii::ExcMessage("Error system::moments::SphericalHarmonic "
                                 "constructor, l_max must be >= 0"));

  for (int group = 0; group < total_groups; ++group) {
    for (int l = 0; l <= max_harmonic_l; ++l) {
      for (int m = -l; m <= l; ++m) {
        moments_.insert(
            std::make_pair<MomentIndex, MomentVector>({group, l, m}, {}));
      }
    }
  }
}
MomentsMap::const_iterator SphericalHarmonic::cbegin() const {
  return moments_.cbegin();
}
MomentsMap::iterator SphericalHarmonic::begin() {
  return moments_.begin();
}
MomentsMap::const_iterator SphericalHarmonic::cend() const {
  return moments_.cend();
}
MomentsMap::iterator SphericalHarmonic::end() {
  return moments_.end();
}
} // namespace moments

} // namespace system

} // namespace bart