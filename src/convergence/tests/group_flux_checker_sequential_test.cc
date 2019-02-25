#include "../group_flux_checker_sequential.h"

#include <memory>

#include <deal.II/lac/petsc_parallel_vector.h>

#include "../../test_helpers/test_helper_functions.h"
#include "../../test_helpers/gmock_wrapper.h"
#include "flux_checker_mock.h"

using ::testing::_;

class GroupFluxCheckerSequentialTest : public ::testing::Test {
 protected:

  bart::convergence::GroupFluxCheckerSequential sequential_tester;
  
  std::unique_ptr<bart::convergence::FluxCheckerI> tester_ptr;
  std::unique_ptr<bart::convergence::FluxCheckerMock> tester_mock;

  bart::data::ScalarGroupFluxPtrs current;
  bart::data::ScalarGroupFluxPtrs previous;

  bart::data::AngularGroupFluxPtrs current_angular;
  bart::data::AngularGroupFluxPtrs previous_angular;

  void SetUp() override;
  void FillGroupFluxes(bart::data::ScalarGroupFluxPtrs &to_fill, int n_groups);
  void FillGroupFluxes(bart::data::AngularGroupFluxPtrs &to_fill,
                       int n_groups, int n_angles);

  void MocksToPointers() {
    tester_ptr = std::move(tester_mock);
  };
};

void GroupFluxCheckerSequentialTest::FillGroupFluxes(
    bart::data::ScalarGroupFluxPtrs &to_fill, int n_groups) {
  for (int i = 0; i < n_groups; ++i) {
    auto flux = std::make_unique<bart::data::Flux>();
    flux->reinit(MPI_COMM_WORLD, 5, 5);
    auto random_vector = btest::RandomVector(5, 0, 2);
    for (unsigned int j = 0; j < flux->size(); ++j)
      (*flux)(j) = random_vector[j];
    flux->compress(dealii::VectorOperation::values::insert);
    to_fill[i] = std::move(flux);
  }
}

void GroupFluxCheckerSequentialTest::FillGroupFluxes(
    bart::data::AngularGroupFluxPtrs &to_fill, int n_groups, int n_angles) {
  for (int i = 0; i < n_groups; ++i) {
    for (int j = 0; j < n_angles; ++j) {
      auto flux = std::make_unique<bart::data::Flux>();
      flux->reinit(MPI_COMM_WORLD, 5, 5);
      auto random_vector = btest::RandomVector(5, 0, 2);
      for (unsigned int j = 0; j < flux->size(); ++j)
        (*flux)(j) = random_vector[j];
      flux->compress(dealii::VectorOperation::values::insert);
      to_fill[std::make_pair(i,j)] = std::move(flux);
    }
  }
}

void GroupFluxCheckerSequentialTest::SetUp() {
  tester_mock = std::make_unique<bart::convergence::FluxCheckerMock>();
}

// Tests where there are no mock calls
class GroupFluxCheckerSeqTestEmptyMock : public GroupFluxCheckerSequentialTest {
 protected:
  bart::convergence::GroupFluxCheckerSequential sequential_tester;
  void SetUp() override {
    MocksToPointers();
    sequential_tester.ProvideChecker(tester_ptr);
  }
};

TEST_F(GroupFluxCheckerSeqTestEmptyMock, Constructor) {
  EXPECT_EQ(tester_ptr, nullptr);
}

TEST_F(GroupFluxCheckerSeqTestEmptyMock, EmptyGroups) {
  EXPECT_ANY_THROW(sequential_tester.CheckIfConverged(current, previous));
}

TEST_F(GroupFluxCheckerSeqTestEmptyMock, DifferentGroupSizes) {
  FillGroupFluxes(current, 5);
  FillGroupFluxes(previous, 4);
  FillGroupFluxes(current_angular, 3, 2);
  FillGroupFluxes(previous_angular, 4, 2);
  EXPECT_ANY_THROW(sequential_tester.CheckIfConverged(current, previous));
  EXPECT_ANY_THROW(sequential_tester.CheckIfConverged(current_angular,
                                                 previous_angular));
}

TEST_F(GroupFluxCheckerSequentialTest, GoodMatch) {
  EXPECT_CALL(*tester_mock, CheckIfConverged(_,_)).
      WillRepeatedly(::testing::Return(true));
  
  MocksToPointers();
  sequential_tester.ProvideChecker(tester_ptr);
  
  FillGroupFluxes(current, 3);
  FillGroupFluxes(previous, 3);

  FillGroupFluxes(current_angular, 3, 2);
  FillGroupFluxes(previous_angular, 3, 2);

  EXPECT_TRUE(sequential_tester.CheckIfConverged(current, previous));
  EXPECT_TRUE(sequential_tester.CheckIfConverged(current_angular,
                                                 previous_angular));
}

TEST_F(GroupFluxCheckerSequentialTest, BadGroupMatch) {
  EXPECT_CALL(*tester_mock, CheckIfConverged(_,_)).
      WillRepeatedly(::testing::Return(true));
  
  MocksToPointers();
  sequential_tester.ProvideChecker(tester_ptr);
  
  FillGroupFluxes(current, 3);
  FillGroupFluxes(previous, 3);

  EXPECT_ANY_THROW(sequential_tester.CheckIfConverged(current, previous));
}

TEST_F(GroupFluxCheckerSequentialTest, BadMatch) {
  EXPECT_CALL(*tester_mock, CheckIfConverged(_,_)).
      WillOnce(::testing::Return(true)).
      WillOnce(::testing::Return(false));
  
  MocksToPointers();
  sequential_tester.ProvideChecker(tester_ptr);
  
  FillGroupFluxes(current, 3);
  FillGroupFluxes(previous, 3);

  EXPECT_FALSE(sequential_tester.CheckIfConverged(current, previous));
  EXPECT_EQ(sequential_tester.GetFailedGroup(), 1);
}

TEST_F(GroupFluxCheckerSequentialTest, BadAngularMatch) {
  EXPECT_CALL(*tester_mock, CheckIfConverged(_,_)).
      WillOnce(::testing::Return(true)).
      WillOnce(::testing::Return(true)).
      WillOnce(::testing::Return(true)).
      WillOnce(::testing::Return(false));
  
  MocksToPointers();
  sequential_tester.ProvideChecker(tester_ptr);
  
  FillGroupFluxes(current_angular, 3, 2);
  FillGroupFluxes(previous_angular, 3, 2);

  EXPECT_FALSE(sequential_tester.CheckIfConverged(current_angular,
                                                  previous_angular));
  EXPECT_EQ(sequential_tester.GetFailedGroup(), 1);
}
