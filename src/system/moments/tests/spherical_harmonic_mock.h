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
  MOCK_METHOD(int, total_groups, (), (const, override));
  MOCK_METHOD(int, max_harmonic_l, (), (const, override));
  MOCK_METHOD(const MomentsMap&, moments, (), (const, override));
  MOCK_METHOD(const MomentVector&, GetMoment, (const MomentIndex), (const, override));
  MOCK_METHOD(MomentsMap::const_iterator, cbegin, (), (const, override));
  MOCK_METHOD(MomentsMap::iterator, begin, (), (override));
  MOCK_METHOD(MomentsMap::const_iterator, cend, (), (const, override));
  MOCK_METHOD(MomentsMap::iterator, end, (), (override));

  MOCK_METHOD(const MomentVector&, BracketOp, (const MomentIndex), (const));
  MOCK_METHOD(MomentVector&, BracketOp, (const MomentIndex));
  
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