#include "system/system_functions.h"

#include "domain/tests/definition_mock.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "system/solution/mpi_group_angular_solution.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_assertions.h"
#include "test_helpers/test_helper_functions.h"

#include "system/terms/term.h"
#include "system/moments/spherical_harmonic.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "system/terms/tests/linear_term_mock.h"
#include "system/terms/tests/bilinear_term_mock.h"
#include "system/solution/solution_types.h"

namespace  {

using namespace bart;
using ::testing::DoDefault;
using ::testing::Return, ::testing::ReturnRef;
using ::testing::WhenDynamicCastTo, ::testing::NotNull;
using ::testing::_, ::testing::NiceMock;

void StampMPIVector(bart::system::MPIVector &to_fill, double value = 2) {
  auto [local_begin, local_end] = to_fill.local_range();
  for (unsigned int i = local_begin; i < local_end; ++i)
    to_fill(i) += value;
  to_fill.compress(dealii::VectorOperation::add);
}

template <typename DimensionWrapper>
class SystemFunctionsSetUpMPIAngularSolutionTests :
    public ::testing::Test,
    public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
  bart::system::solution::MPIGroupAngularSolutionMock mock_solution;
  domain::DefinitionMock<dim> mock_definition;

  const int total_angles_ = 3;
  std::map<bart::system::AngleIndex, bart::system::MPIVector> solution_map_;
  void SetUp() override;
};

template <typename DimensionWrapper>
void SystemFunctionsSetUpMPIAngularSolutionTests<DimensionWrapper>::SetUp() {
  this->SetUpDealii();
  ON_CALL(mock_solution, total_angles()).WillByDefault(Return(total_angles_));

  for (int i = 0; i < total_angles_; ++i) {
    bart::system::MPIVector mpi_vector;
    solution_map_.insert_or_assign(i, mpi_vector);
  }
  ON_CALL(mock_solution, solutions()).WillByDefault(ReturnRef(solution_map_));
  ON_CALL(mock_definition, locally_owned_dofs())
      .WillByDefault(Return(this->locally_owned_dofs_));
}

TYPED_TEST_SUITE(SystemFunctionsSetUpMPIAngularSolutionTests,
    bart::testing::AllDimensions);

TYPED_TEST(SystemFunctionsSetUpMPIAngularSolutionTests, BadNangles) {
  constexpr int dim = this->dim;

  std::array<int, 4> bad_total_angles{0, -1, 2, 4};

  for (const auto angle : bad_total_angles) {
    EXPECT_CALL(this->mock_solution, total_angles()).WillOnce(Return(angle));
    EXPECT_CALL(this->mock_solution, solutions()).WillOnce(DoDefault());
    EXPECT_ANY_THROW({
      bart::system::SetUpMPIAngularSolution<dim>(this->mock_solution,
                                                 this->mock_definition);
                     });
  }
}

TYPED_TEST(SystemFunctionsSetUpMPIAngularSolutionTests, SetUpDefaultValue) {
  constexpr int dim = this->dim;
  EXPECT_CALL(this->mock_solution, total_angles()).WillOnce(DoDefault());
  EXPECT_CALL(this->mock_solution, solutions()).WillOnce(DoDefault());
  EXPECT_CALL(this->mock_definition, locally_owned_dofs()).WillOnce(DoDefault());
  EXPECT_NO_THROW({
    bart::system::SetUpMPIAngularSolution<dim>(this->mock_solution,
                                               this->mock_definition);
  });
  bart::system::MPIVector expected_vector;
  expected_vector.reinit(this->locally_owned_dofs_, MPI_COMM_WORLD);
  StampMPIVector(expected_vector, 1.0);

  for (const auto& solution : this->solution_map_) {
    auto& mpi_vector = solution.second;
    ASSERT_GT(mpi_vector.size(), 0);
    EXPECT_TRUE(bart::test_helpers::CompareMPIVectors(expected_vector, mpi_vector));
  }
}

TYPED_TEST(SystemFunctionsSetUpMPIAngularSolutionTests, ProvidedValues) {
  constexpr int dim = this->dim;
  const double value_to_set = test_helpers::RandomDouble(0, 20);

  EXPECT_CALL(this->mock_solution, total_angles()).WillOnce(DoDefault());
  EXPECT_CALL(this->mock_solution, solutions()).WillOnce(DoDefault());
  EXPECT_CALL(this->mock_definition, locally_owned_dofs()).WillOnce(DoDefault());
  EXPECT_NO_THROW({
    bart::system::SetUpMPIAngularSolution<dim>(this->mock_solution,
                                               this->mock_definition,
                                               value_to_set);
                   });
  bart::system::MPIVector expected_vector;
  expected_vector.reinit(this->locally_owned_dofs_, MPI_COMM_WORLD);
  expected_vector = 0;
  StampMPIVector(expected_vector, value_to_set);

  for (const auto& solution : this->solution_map_) {
    auto& mpi_vector = solution.second;
    ASSERT_GT(mpi_vector.size(), 0);
    EXPECT_TRUE(bart::test_helpers::CompareMPIVectors(expected_vector, mpi_vector));
  }
}

// == InitializeSystem Tests ===================================================

class SystemFunctionsInitializeSystemTest : public ::testing::Test {
 public:
  bart::system::System test_system;
};

TEST_F(SystemFunctionsInitializeSystemTest, DefaultCall) {
  using VariableLinearTerms = system::terms::VariableLinearTerms;
  using ExpectedRHSType = bart::system::terms::MPILinearTerm;
  using ExpectedLHSType = bart::system::terms::MPIBilinearTerm;
  using ExpectedMomentsType = bart::system::moments::SphericalHarmonic;

  const int total_groups = bart::test_helpers::RandomDouble(1, 10);
  const int total_angles = total_groups + 1;
  std::unordered_set<VariableLinearTerms> source_terms{
      VariableLinearTerms::kScatteringSource,
      VariableLinearTerms::kFissionSource};

  system::InitializeSystem(test_system, total_groups, total_angles);

  EXPECT_EQ(test_system.total_angles, total_angles);
  EXPECT_EQ(test_system.total_groups, total_groups);
  EXPECT_EQ(test_system.k_effective, 1.0);
  ASSERT_THAT(test_system.right_hand_side_ptr_.get(),
              WhenDynamicCastTo<ExpectedRHSType *>(NotNull()));
  ASSERT_THAT(test_system.left_hand_side_ptr_.get(),
              WhenDynamicCastTo<ExpectedLHSType *>(NotNull()));
  EXPECT_EQ(test_system.right_hand_side_ptr_->GetVariableTerms(), source_terms);
  ASSERT_THAT(test_system.current_moments.get(),
              WhenDynamicCastTo<ExpectedMomentsType *>(NotNull()));
  ASSERT_THAT(test_system.previous_moments.get(),
              WhenDynamicCastTo<ExpectedMomentsType *>(NotNull()));
  EXPECT_EQ(test_system.current_moments->moments().size(), total_groups);
  EXPECT_EQ(test_system.previous_moments->moments().size(), total_groups);
}

TEST_F(SystemFunctionsInitializeSystemTest, NonEigenvalueProblem) {
  using VariableLinearTerms = system::terms::VariableLinearTerms;
  using ExpectedRHSType = bart::system::terms::MPILinearTerm;
  using ExpectedLHSType = bart::system::terms::MPIBilinearTerm;
  using ExpectedMomentsType = bart::system::moments::SphericalHarmonic;

  const int total_groups = bart::test_helpers::RandomDouble(1, 10);
  const int total_angles = total_groups + 1;
  std::unordered_set<VariableLinearTerms> source_terms{
      VariableLinearTerms::kScatteringSource};

  system::InitializeSystem(test_system, total_groups, total_angles, false);

  EXPECT_EQ(test_system.total_angles, total_angles);
  EXPECT_EQ(test_system.total_groups, total_groups);
  EXPECT_EQ(test_system.k_effective, std::nullopt);
  ASSERT_THAT(test_system.right_hand_side_ptr_.get(),
              WhenDynamicCastTo<ExpectedRHSType *>(NotNull()));
  ASSERT_THAT(test_system.left_hand_side_ptr_.get(),
              WhenDynamicCastTo<ExpectedLHSType *>(NotNull()));
  EXPECT_EQ(test_system.right_hand_side_ptr_->GetVariableTerms(), source_terms);
  ASSERT_THAT(test_system.current_moments.get(),
              WhenDynamicCastTo<ExpectedMomentsType *>(NotNull()));
  ASSERT_THAT(test_system.previous_moments.get(),
              WhenDynamicCastTo<ExpectedMomentsType *>(NotNull()));
  EXPECT_EQ(test_system.current_moments->moments().size(), total_groups);
  EXPECT_EQ(test_system.previous_moments->moments().size(), total_groups);
}

TEST_F(SystemFunctionsInitializeSystemTest, BoundaryConditions) {
  using VariableLinearTerms = system::terms::VariableLinearTerms;
  using ExpectedRHSType = bart::system::terms::MPILinearTerm;
  using ExpectedLHSType = bart::system::terms::MPIBilinearTerm;
  using ExpectedMomentsType = bart::system::moments::SphericalHarmonic;

  const int total_groups = bart::test_helpers::RandomDouble(1, 10);
  const int total_angles = total_groups + 1;
  std::unordered_set<VariableLinearTerms> source_terms{
      VariableLinearTerms::kScatteringSource,
      VariableLinearTerms::kReflectiveBoundaryCondition};

  system::InitializeSystem(test_system, total_groups, total_angles, false, true);

  EXPECT_EQ(test_system.total_angles, total_angles);
  EXPECT_EQ(test_system.total_groups, total_groups);
  EXPECT_EQ(test_system.k_effective, std::nullopt);
  ASSERT_THAT(test_system.right_hand_side_ptr_.get(),
              WhenDynamicCastTo<ExpectedRHSType *>(NotNull()));
  ASSERT_THAT(test_system.left_hand_side_ptr_.get(),
              WhenDynamicCastTo<ExpectedLHSType *>(NotNull()));
  EXPECT_EQ(test_system.right_hand_side_ptr_->GetVariableTerms(), source_terms);
  ASSERT_THAT(test_system.current_moments.get(),
              WhenDynamicCastTo<ExpectedMomentsType *>(NotNull()));
  ASSERT_THAT(test_system.previous_moments.get(),
              WhenDynamicCastTo<ExpectedMomentsType *>(NotNull()));
  EXPECT_EQ(test_system.current_moments->moments().size(), total_groups);
  EXPECT_EQ(test_system.previous_moments->moments().size(), total_groups);
}

TEST_F(SystemFunctionsInitializeSystemTest, ErrorOnSecondCall) {
  using VariableLinearTerms = system::terms::VariableLinearTerms;
  using ExpectedRHSType = bart::system::terms::MPILinearTerm;
  using ExpectedLHSType = bart::system::terms::MPIBilinearTerm;
  using ExpectedMomentsType = bart::system::moments::SphericalHarmonic;

  const int total_groups = bart::test_helpers::RandomDouble(1, 10);
  const int total_angles = total_groups + 1;
  std::unordered_set<VariableLinearTerms> source_terms{
      VariableLinearTerms::kScatteringSource,
      VariableLinearTerms::kFissionSource};

  system::InitializeSystem(test_system, total_groups, total_angles);
  EXPECT_ANY_THROW(bart::system::InitializeSystem(test_system, total_groups,
                                                  total_angles));
}

// == SetUpSystemTerms Tests ===================================================

template <typename DimensionWrapper>
class SystemFunctionsSetUpSystemTermsTests : public ::testing::Test,
                                             bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using DomainType = domain::DefinitionMock<dim>;
  using RhsTermType = bart::system::terms::LinearTermMock;
  using LhsTermType = bart::system::terms::BilinearTermMock;
  using VariableLinearTerms = bart::system::terms::VariableLinearTerms;

  bart::system::System test_system;
  std::shared_ptr<domain::DefinitionI<dim>> definition_ptr;
  DomainType* domain_mock_obs_ptr_;
  RhsTermType* rhs_mock_obs_ptr_;
  LhsTermType* lhs_mock_obs_ptr_;

  std::shared_ptr<bart::system::MPISparseMatrix> system_matrix_ptr_;
  std::shared_ptr<bart::system::MPIVector> system_vector_ptr_;

  std::unordered_set<VariableLinearTerms> source_terms_{
      VariableLinearTerms::kScatteringSource,
      VariableLinearTerms::kFissionSource};

  void SetUp() override;
};

template <typename DimensionWrapper>
void SystemFunctionsSetUpSystemTermsTests<DimensionWrapper>::SetUp() {
  this->SetUpDealii();

  test_system.total_groups = bart::test_helpers::RandomDouble(2, 4);
  test_system.total_angles = test_system.total_groups - 1;
  test_system.right_hand_side_ptr_ = std::make_unique<RhsTermType>();
  test_system.left_hand_side_ptr_ = std::make_unique<LhsTermType>();

  rhs_mock_obs_ptr_ = dynamic_cast<RhsTermType*>(test_system.right_hand_side_ptr_.get());
  lhs_mock_obs_ptr_ = dynamic_cast<LhsTermType*>(test_system.left_hand_side_ptr_.get());

  ON_CALL(*rhs_mock_obs_ptr_, GetVariableTerms())
      .WillByDefault(Return(source_terms_));

  definition_ptr = std::make_shared<DomainType>();
  domain_mock_obs_ptr_ = dynamic_cast<DomainType*>(definition_ptr.get());

  system_matrix_ptr_ = std::make_shared<bart::system::MPISparseMatrix>();
  system_matrix_ptr_->reinit(this->matrix_1);
  system_vector_ptr_ = std::make_shared<bart::system::MPIVector>();
  system_vector_ptr_->reinit(this->vector_1);

  ON_CALL(*domain_mock_obs_ptr_, MakeSystemMatrix())
      .WillByDefault(Return(system_matrix_ptr_));
  ON_CALL(*domain_mock_obs_ptr_, MakeSystemVector())
      .WillByDefault(Return(system_vector_ptr_));
}

TYPED_TEST_SUITE(SystemFunctionsSetUpSystemTermsTests, bart::testing::AllDimensions);

TYPED_TEST(SystemFunctionsSetUpSystemTermsTests, SetUpProperly) {
  auto& test_system = this->test_system;
  const int total_groups = test_system.total_groups;
  const int total_angles = test_system.total_angles;

  EXPECT_CALL(*this->rhs_mock_obs_ptr_, GetVariableTerms())
      .WillOnce(DoDefault());

  EXPECT_CALL(*this->domain_mock_obs_ptr_, MakeSystemMatrix())
      .Times(total_angles * total_groups)
      .WillRepeatedly(DoDefault());
  EXPECT_CALL(*this->domain_mock_obs_ptr_, MakeSystemVector())
      .Times(total_groups * total_angles * (1 + this->source_terms_.size()))
      .WillRepeatedly(DoDefault());

  for (int group = 0; group < total_groups; ++group) {
    for (int angle = 0; angle < total_angles; ++angle) {
      bart::system::Index index{group, angle};
      EXPECT_CALL(*this->rhs_mock_obs_ptr_, SetFixedTermPtr(index, NotNull()));
      EXPECT_CALL(*this->lhs_mock_obs_ptr_, SetFixedTermPtr(index, NotNull()));
      for (auto term : this->source_terms_)
        EXPECT_CALL(*this->rhs_mock_obs_ptr_, SetVariableTermPtr(index, term, NotNull()));
    }
  }

  bart::system::SetUpSystemTerms(test_system, *this->definition_ptr);
}

// ===== SetUpSystemMomentsTests ===============================================

class SystemFunctionsSetUpSystemMomentsTests : public ::testing::Test {
 public:
  using MomentsType = NiceMock<bart::system::moments::SphericalHarmonicMock>;

  bart::system::System test_system;
  MomentsType* current_moments_obs_ptr_;
  MomentsType* previous_moments_obs_ptr_;

  template <typename T> inline MomentsType* MockCast(T* to_cast) {
    return dynamic_cast<MomentsType*>(to_cast); }

  const int solution_size = test_helpers::RandomDouble(1, 100);
  bart::system::moments::MomentsMap current_moments_, previous_moments_;

  void SetUp() override;
};

void SystemFunctionsSetUpSystemMomentsTests::SetUp() {
  test_system.current_moments = std::make_unique<MomentsType>();
  test_system.previous_moments = std::make_unique<MomentsType>();

  current_moments_obs_ptr_ = MockCast(test_system.current_moments.get());
  previous_moments_obs_ptr_ = MockCast(test_system.previous_moments.get());

  const int n_groups = bart::test_helpers::RandomDouble(1, 4);
  const int max_harmonic_l = bart::test_helpers::RandomDouble(0, 3);
  test_system.total_groups = n_groups;

  for (auto& mock_moment_ptr : {current_moments_obs_ptr_,
                                previous_moments_obs_ptr_}) {
    ON_CALL(*mock_moment_ptr, total_groups())
        .WillByDefault(Return(n_groups));
    ON_CALL(*mock_moment_ptr, max_harmonic_l())
        .WillByDefault(Return(max_harmonic_l));
  }


  for (int group = 0; group < n_groups; ++group) {
    for (int harmonic_l = 0; harmonic_l <= max_harmonic_l; ++harmonic_l) {
      for (int harmonic_m = -harmonic_l; harmonic_m <= harmonic_l; ++harmonic_m) {
        system::moments::MomentIndex index{group, harmonic_l, harmonic_m};
        current_moments_.emplace(index, dealii::Vector<double>{});
        previous_moments_.emplace(index, dealii::Vector<double>{});
      }
    }
  }

  ON_CALL(*current_moments_obs_ptr_, moments())
      .WillByDefault(ReturnRef(current_moments_));
  ON_CALL(*current_moments_obs_ptr_, begin())
      .WillByDefault(Return(current_moments_.begin()));
  ON_CALL(*current_moments_obs_ptr_, end())
      .WillByDefault(Return(current_moments_.end()));
  ON_CALL(*previous_moments_obs_ptr_, moments())
      .WillByDefault(ReturnRef(previous_moments_));
  ON_CALL(*previous_moments_obs_ptr_, begin())
      .WillByDefault(Return(previous_moments_.begin()));
  ON_CALL(*previous_moments_obs_ptr_, end())
      .WillByDefault(Return(previous_moments_.end()));
}

TEST_F(SystemFunctionsSetUpSystemMomentsTests, SetUpProperly) {
  for (auto& mock_obs_ptr : {current_moments_obs_ptr_,
                             previous_moments_obs_ptr_}) {
    EXPECT_CALL(*mock_obs_ptr, begin()).WillOnce(DoDefault());
    EXPECT_CALL(*mock_obs_ptr, end()).WillOnce(DoDefault());
  }

  system::SetUpSystemMoments(test_system, solution_size);
  dealii::Vector<double> expected(solution_size);
  expected = 1;

  for (const auto& moment_map : {current_moments_, previous_moments_}) {
    for (const auto& moment_pair : moment_map) {
      const auto& moment_vector = moment_pair.second;
      ASSERT_EQ(moment_vector.size(), solution_size);
      EXPECT_EQ(moment_vector, expected);
    }
  }
}

// ===== SetUpSystemAngularSolution Tests ======================================

class SystemFunctionsSetUpEnergyGroupToAngularSolutionPtrMapIntTests
    : public ::testing::Test {
 public:
  system::solution::EnergyGroupToAngularSolutionPtrMap solution_map_;

  const int total_groups_ = test_helpers::RandomDouble(2, 4);
  const int total_angles_{total_groups_ + 1};

};

TEST_F(SystemFunctionsSetUpEnergyGroupToAngularSolutionPtrMapIntTests,
       DefaultCall) {
  system::SetUpEnergyGroupToAngularSolutionPtrMap(solution_map_,
                                                  total_groups_,
                                                  total_angles_);

  using ExpectedType = bart::system::solution::MPIGroupAngularSolution;
  std::vector<int> groups{};
  EXPECT_EQ(solution_map_.size(), total_groups_);
  for (auto& [energy_group, solution_ptr] : solution_map_) {
    ASSERT_THAT(solution_ptr.get(),
                WhenDynamicCastTo<ExpectedType*>(NotNull()));
    EXPECT_EQ(solution_ptr->total_angles(), total_angles_);
    EXPECT_LT(energy_group.get(), total_groups_);
    EXPECT_GE(energy_group.get(), 0);
    EXPECT_EQ(std::count(groups.cbegin(), groups.cend(), energy_group.get()), 0);
    groups.push_back(energy_group.get());
  }
}



} // namespace