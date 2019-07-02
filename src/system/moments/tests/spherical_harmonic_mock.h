#ifndef BART_SRC_SYSTEM_MOMENTS_TESTS_SPHERICAL_HARMONIC_MOCK_H_
#define BART_SRC_SYSTEM_MOMENTS_TESTS_SPHERICAL_HARMONIC_MOCK_H_

#include "system/moments/spherical_harmonic_i.h"
#include "system/moments/spherical_harmonic_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace system {

namespace moments {

class SphericalHarmonicMock : public SphericalHarmonicI {
 public:
  MOCK_CONST_METHOD0(total_groups, int());
  MOCK_CONST_METHOD0(max_harmonic_l, int());
  MOCK_CONST_METHOD0(moments, const MomentsMap&());
  MOCK_CONST_METHOD1(GetMoment, const MomentVector&(const MomentIndex));
  MOCK_CONST_METHOD1(BracketOp, const MomentVector&(const MomentIndex));
  MOCK_METHOD1(BracketOp, MomentVector&(const MomentIndex));

  const MomentVector& operator[](const MomentIndex index) const override {
    return BracketOp(index);
  };

  MomentVector& operator[](const MomentIndex index) override {
    return BracketOp(index);
  };
};

} // namespace moments

} // namespace system


} // namespace bart

#endif // BART_SRC_SYSTEM_MOMENTS_TESTS_SPHERICAL_HARMONIC_MOCK_H_