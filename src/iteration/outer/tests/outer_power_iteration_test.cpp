#include "iteration/outer/outer_power_iteration.hpp"

#include <memory>

#include "instrumentation/tests/instrument_mock.h"
#include "iteration/group/tests/group_solve_iteration_mock.hpp"
#include "iteration/subroutine/tests/subroutine_mock.hpp"
#include "eigenvalue/k_eigenvalue/tests/k_eigenvalue_calculator_mock.hpp"
#include "convergence/tests/iteration_completion_checker_mock.hpp"
#include "formulation/updater/tests/fission_source_updater_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "system/system.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "system/terms/tests/linear_term_mock.hpp"

namespace  {

using namespace bart;

using ::testing::A, ::testing::AtLeast, ::testing::Expectation;
using ::testing::Ref, ::testing::Return, ::testing::Sequence, ::testing::_;

/* This fixture tests the operation of the OuterPowerIteration class. This is a mediator class so the tests will verify
 * proper mediation between the dependencies and exposure of data to instruments.
 *
 * At the completion of SetUp() the test_iterator pointee object is set up with mock dependencies accessible via
 * the provided observation pointers. Mock instrumentation is accessible via shared pointers.
 * */
class IterationOuterPowerIterationTest : public ::testing::Test {
 protected:
  using GroupIterator = iteration::group::GroupSolveIterationMock;
  using ConvergenceChecker = convergence::IterationCompletionCheckerMock<double>;
  using ConvergenceInstrumentType = instrumentation::InstrumentMock<convergence::Status>;
  using ErrorInstrumentType = instrumentation::InstrumentMock<std::pair<int, double>>;
  using FissionSourceInstrumentMock = instrumentation::InstrumentMock<std::unordered_map<int, dealii::Vector<double>>>;
  using K_EffectiveUpdater = eigenvalue::k_eigenvalue::K_EigenvalueCalculatorMock;
  using OuterPowerIteration = iteration::outer::OuterPowerIteration;
  using MomentsMock = system::moments::SphericalHarmonicMock;
  using RightHandSideMock = system::terms::LinearTermMock;
  using ScalarFluxInstrumentMock = instrumentation::InstrumentMock<std::unordered_map<int, dealii::Vector<double>>>;
  using ScatteringSourceInstrumentMock = instrumentation::InstrumentMock<dealii::Vector<double>>;
  using SourceUpdater = formulation::updater::FissionSourceUpdaterMock;
  using StatusInstrumentType = instrumentation::InstrumentMock<std::string>;
  using SolutionMomentsInstrumentMock = instrumentation::InstrumentMock<bart::system::moments::SphericalHarmonicI>;
  using Subroutine = iteration::subroutine::SubroutineMock;

  std::unique_ptr<OuterPowerIteration> test_iterator;

  // Dependencies
  std::shared_ptr<SourceUpdater> source_updater_ptr_{ std::make_shared<SourceUpdater>() };

  // Mock instruments
  std::shared_ptr<ConvergenceInstrumentType> convergence_instrument_ptr_{ std::make_shared<ConvergenceInstrumentType>() };
  std::shared_ptr<ErrorInstrumentType> error_instrument_ptr_{ std::make_shared<ErrorInstrumentType>() };
  std::shared_ptr<FissionSourceInstrumentMock> fission_source_instrument_ptr_{
      std::make_shared<FissionSourceInstrumentMock>() };
  std::shared_ptr<StatusInstrumentType> status_instrument_ptr_{ std::make_shared<StatusInstrumentType>() };
  std::shared_ptr<ScalarFluxInstrumentMock> scalar_flux_instrument_ptr_{ std::make_shared<ScalarFluxInstrumentMock>() };
  std::shared_ptr<ScatteringSourceInstrumentMock> scattering_source_instrument_ptr_{
      std::make_shared<ScatteringSourceInstrumentMock>() };
  std::shared_ptr<SolutionMomentsInstrumentMock> solution_moments_instrument_ptr_{
    std::make_shared<SolutionMomentsInstrumentMock>() };

  // Supporting objects
  system::System test_system;

  // Mocks observation pointers
  std::shared_ptr<MomentsMock> current_moments_mock_ptr_{std::make_shared<MomentsMock>() };
  GroupIterator* group_iterator_obs_ptr_;
  ConvergenceChecker* convergence_checker_obs_ptr_;
  K_EffectiveUpdater* k_effective_updater_obs_ptr_;
  RightHandSideMock* right_hand_side_obs_ptr_;
  Subroutine* post_iteration_subroutine_obs_ptr_;

  // Test parameters
  static constexpr int total_groups{ 2 };
  static constexpr int total_angles{ 3 };
  static constexpr int iterations_{ 4 };
  static constexpr int total_degrees_of_freedom_{ 10 };

  void SetUp() override;
};

void IterationOuterPowerIterationTest::SetUp() {
  // Dependencies
  auto group_iterator_ptr = std::make_unique<GroupIterator>();
  group_iterator_obs_ptr_ = group_iterator_ptr.get();
  auto convergenge_checker_ptr = std::make_unique<ConvergenceChecker>();
  convergence_checker_obs_ptr_ = convergenge_checker_ptr.get();
  auto k_effective_updater_ptr = std::make_unique<K_EffectiveUpdater>();
  k_effective_updater_obs_ptr_ = k_effective_updater_ptr.get();
  auto post_iteration_subroutine_ptr = std::make_unique<Subroutine>();
  post_iteration_subroutine_obs_ptr_ = post_iteration_subroutine_ptr.get();
  auto right_hand_side_ptr = std::make_unique<RightHandSideMock>();
  right_hand_side_obs_ptr_ = right_hand_side_ptr.get();

  // Set up system
  test_system.total_angles = total_angles;
  test_system.total_groups = total_groups;
  test_system.current_moments = this->current_moments_mock_ptr_;
  test_system.right_hand_side_ptr_ = std::move(right_hand_side_ptr);

  // Construct test object
  test_iterator = std::make_unique<OuterPowerIteration>(
      std::move(group_iterator_ptr),
      std::move(convergenge_checker_ptr),
      std::move(k_effective_updater_ptr),
      source_updater_ptr_);
  test_iterator->AddPostIterationSubroutine(std::move(post_iteration_subroutine_ptr));

  using ConvergenceStatusPort = iteration::outer::data_names::ConvergenceStatusPort;
  instrumentation::GetPort<ConvergenceStatusPort>(*test_iterator).AddInstrument(convergence_instrument_ptr_);
  using FissionSourcePort = iteration::outer::data_names::FissionSourcePort;
  instrumentation::GetPort<FissionSourcePort>(*test_iterator).AddInstrument(fission_source_instrument_ptr_);
  using StatusPort = iteration::outer::data_names::StatusPort;
  instrumentation::GetPort<StatusPort>(*test_iterator).AddInstrument(status_instrument_ptr_);
  using IterationErrorPort = iteration::outer::data_names::IterationErrorPort;
  instrumentation::GetPort<IterationErrorPort>(*test_iterator).AddInstrument(error_instrument_ptr_);
  using ScalarFluxPort = iteration::outer::data_names::ScalarFluxPort;
  instrumentation::GetPort<ScalarFluxPort>(*test_iterator).AddInstrument(scalar_flux_instrument_ptr_);
  using ScatteringSourcePort = iteration::outer::data_names::ScatteringSourcePort;
  instrumentation::GetPort<ScatteringSourcePort>(*test_iterator).AddInstrument(scattering_source_instrument_ptr_);
  using SolutionMomentsPort = iteration::outer::data_names::SolutionMomentsPort;
  instrumentation::GetPort<SolutionMomentsPort>(*test_iterator).AddInstrument(solution_moments_instrument_ptr_);
}

/* Constructor (called in SetUp()) should have stored the correct pointers in the test object */
TEST_F(IterationOuterPowerIterationTest, Constructor) {
  EXPECT_NE(this->test_iterator, nullptr);
  EXPECT_NE(this->test_iterator->group_iterator_ptr(), nullptr);
  EXPECT_NE(this->test_iterator->source_updater_ptr(), nullptr);
  EXPECT_NE(this->test_iterator->convergence_checker_ptr(), nullptr);
  EXPECT_NE(this->test_iterator->k_effective_updater_ptr(), nullptr);
  EXPECT_NE(this->test_iterator->post_iteration_subroutine_ptr(), nullptr);
  EXPECT_EQ(this->source_updater_ptr_.use_count(), 2);
}

/* Constructor should throw an error if pointers passed to dependencies are null. */
TEST_F(IterationOuterPowerIterationTest, ConstructorErrors) {

  for (int i = 0; i < 4; ++i) {
    auto convergence_checker_ptr = (i == 0) ? nullptr :
        std::make_unique<convergence::IterationCompletionCheckerMock<double>>();
    auto k_effective_updater_ptr = (i == 1) ? nullptr :
        std::make_unique<eigenvalue::k_eigenvalue::K_EigenvalueCalculatorMock>();
    auto source_updater_ptr = (i == 2) ? nullptr : this->source_updater_ptr_;
    auto group_iterator_ptr = (i == 3) ? nullptr :
        std::make_unique<iteration::group::GroupSolveIterationMock>();

    EXPECT_ANY_THROW({
                       iteration::outer::OuterPowerIteration test_iterator(
                           std::move(group_iterator_ptr),
                           std::move(convergence_checker_ptr),
                           std::move(k_effective_updater_ptr),
                           source_updater_ptr);
                     });
  }
}

/* A call to Iterate() should properly mediate the calls to the owned classes, including the group iteration class, and
 * convergence-checker. */
TEST_F(IterationOuterPowerIterationTest, IterateToConvergenceTest) {

  for (int group = 0; group < this->total_groups; ++group) {
    for (int angle = 0; angle < this->total_angles; ++angle) {
      EXPECT_CALL(*this->source_updater_ptr_, UpdateFissionSource(Ref(this->test_system),
                                                                  bart::system::EnergyGroup(group),
                                                                  quadrature::QuadraturePointIndex(angle)))
          .Times(this->iterations_);
    }
  }

  // K_Effective updater return values
  std::array<double, iterations_ + 1> k_effective_by_iteration;
  k_effective_by_iteration.fill(0);
  Sequence k_effective_calls;
  std::vector<double> expected_errors;

  for (int i = 0; i < this->iterations_; ++i) {
    const double iteration_k_effective{i * 1.5 };
    k_effective_by_iteration.at(i + 1) = iteration_k_effective;

    EXPECT_CALL(*this->k_effective_updater_obs_ptr_, CalculateK_Eigenvalue(Ref(this->test_system)))
        .InSequence(k_effective_calls)
        .WillOnce(Return(iteration_k_effective));

    convergence::Status convergence_status;
    convergence_status.is_complete = (i == (this->iterations_ - 1));
    if (i > 0) {
      double expected_error = (iteration_k_effective - k_effective_by_iteration.at(i)) / (iteration_k_effective);
      expected_errors.push_back(expected_error);
      convergence_status.delta = expected_error;
    } else {
      convergence_status.delta = std::nullopt;
    }

    EXPECT_CALL(*this->convergence_checker_obs_ptr_,
                ConvergenceStatus(k_effective_by_iteration.at(i + 1), k_effective_by_iteration.at(i)))
        .WillOnce(Return(convergence_status));
  }

  EXPECT_CALL(*this->convergence_checker_obs_ptr_, Reset());
  EXPECT_CALL(*this->group_iterator_obs_ptr_, Iterate(Ref(this->test_system))).Times(this->iterations_);
  EXPECT_CALL(*this->convergence_instrument_ptr_, Read(_)).Times(this->iterations_);
  EXPECT_CALL(*this->post_iteration_subroutine_obs_ptr_, Execute(Ref(this->test_system))).Times(this->iterations_);
  EXPECT_CALL(*this->status_instrument_ptr_, Read(_)).Times(AtLeast(this->iterations_));
  EXPECT_CALL(*this->error_instrument_ptr_, Read(_)).Times(this->iterations_ - 1);

  using VariableLinearTerms = system::terms::VariableLinearTerms;
  EXPECT_CALL(*this->right_hand_side_obs_ptr_, GetVariableTerms())
      .Times(iterations_)
      .WillRepeatedly(Return(std::unordered_set{VariableLinearTerms::kScatteringSource,
                                                VariableLinearTerms::kFissionSource}));

  // Instrumentation
  std::unordered_map<int, dealii::Vector<double>> flux_map;
  std::unordered_map<int, dealii::Vector<double>> fission_source_map;

  for (int group = 0; group < total_groups; ++group) {
    auto fission_source_vector_ptr = std::make_shared<system::MPIVector>(MPI_COMM_WORLD, total_degrees_of_freedom_,
                                                                         total_degrees_of_freedom_);
    fission_source_vector_ptr->add(2 * (group + 1 ));
    fission_source_vector_ptr->compress(dealii::VectorOperation::add);
    flux_map[group] = dealii::Vector<double>(total_degrees_of_freedom_);
    fission_source_map[group] = dealii::Vector<double>(total_degrees_of_freedom_);
    flux_map.at(group).add(group + 1);
    fission_source_map.at(group).add(2 * (group + 1 ));

    EXPECT_CALL(*this->current_moments_mock_ptr_, GetMoment(std::array{group, 0, 0}))
        .Times(this->iterations_)
        .WillRepeatedly(::testing::ReturnRef(flux_map.at(group)));
    EXPECT_CALL(*this->right_hand_side_obs_ptr_, GetVariableTermPtr(group, VariableLinearTerms::kFissionSource))
        .Times(this->iterations_)
        .WillRepeatedly(Return(fission_source_vector_ptr));
  }
  EXPECT_CALL(*this->scalar_flux_instrument_ptr_, Read(flux_map)).Times(this->iterations_);
  EXPECT_CALL(*this->fission_source_instrument_ptr_, Read(fission_source_map)).Times(iterations_);

  EXPECT_CALL(*this->solution_moments_instrument_ptr_, Read(Ref(*this->current_moments_mock_ptr_)))
      .Times(this->iterations_);

  auto scattering_source_vector_ptr = std::make_shared<system::MPIVector>();
  EXPECT_CALL(*this->right_hand_side_obs_ptr_, GetVariableTermPtr(0, VariableLinearTerms::kScatteringSource))
      .Times(this->iterations_)
      .WillRepeatedly(Return(scattering_source_vector_ptr));

  EXPECT_CALL(*this->scattering_source_instrument_ptr_, Read(_)).Times(iterations_);

  this->test_iterator->IterateToConvergence(this->test_system);
}



} // namespace