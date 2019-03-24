#ifndef BART_SRC_CONVERGENCE_TESTS_SINGLE_CHECKER_TEST_H_
#define BART_SRC_CONVERGENCE_TESTS_SINGLE_CHECKER_TEST_H_

#include "convergence/single_checker.h"

#include <optional>

#include <gtest/gtest.h>

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace testing {

template<typename CompareType>
class SingleCheckerTest : public ::testing::Test {
 protected:
  using SingleChecker = bart::convergence::SingleChecker<CompareType>;

  void TestBaseMethods(SingleChecker *checker) {
    auto delta = std::make_optional<double>(0.185);

    checker->SetMaxDelta(0.185);
    EXPECT_EQ(checker->delta(), delta);
    EXPECT_ANY_THROW(checker->SetMaxDelta(-0.185));
  }

};

} // namespace testing

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_TESTS_SINGLE_CHECKER_TEST_H_