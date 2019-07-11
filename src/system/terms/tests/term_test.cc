#include "system/terms/term.h"

#include "test_helpers/test_assertions.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

/* This testing suite is designed to test the linear and bilinear term classes
 * to ensure they are adding vectors and matrices properly when the
 * GetFullTerm method is called. For tests that verify that the setters and
 * getters work, see linear_term_test.cc. These tests will use a full dealii
 * test environment to provide the needed MPI matrices and vectors.
 */

class SystemTermsFullTermTest : public ::testing::Test,
                                public bart::testing::DealiiTestDomain<2> {
 protected:
  void SetUp() override;
};

void SystemTermsFullTermTest::SetUp() {
  SetUpDealii();
}

TEST_F(SystemTermsFullTermTest, BilinearFullTermOperationMPI) {
  auto other_source = system::terms::VariableBilinearTerms::kOther;
  system::terms::MPIBilinearTerm test_bilinear_term({other_source});

  auto fixed_term_ptr = std::make_shared<system::MPISparseMatrix>();
  auto variable_term_ptr = std::make_shared<system::MPISparseMatrix>();

  fixed_term_ptr->reinit(matrix_1);
  variable_term_ptr->reinit(matrix_2);

  StampMatrix(*fixed_term_ptr, 2);
  StampMatrix(*variable_term_ptr, 1);
  StampMatrix(matrix_3, 3);

  test_bilinear_term.SetFixedTermPtr({0, 0}, fixed_term_ptr);
  test_bilinear_term.SetVariableTermPtr({0, 0}, other_source, variable_term_ptr);

  auto term_matrix_ptr = test_bilinear_term.GetFullTermPtr({0,0});
  EXPECT_TRUE(bart::testing::CompareMPIMatrices(matrix_3, *term_matrix_ptr));
}

TEST_F(SystemTermsFullTermTest, LinearFullTermOperationMPI) {
  using VariableTerms = system::terms::VariableLinearTerms;
  system::terms::MPILinearTerm test_linear_term({VariableTerms::kFissionSource,
                                                 VariableTerms::kOther});

  auto fixed_term_ptr = std::make_shared<system::MPIVector>();
  auto fission_term_ptr = std::make_shared<system::MPIVector>();
  auto other_term_ptr = std::make_shared<system::MPIVector>();

  auto set_value = [&](system::MPIVector& vector, int value) {
    vector.reinit(vector_1);
    vector.add(value);
    vector.compress(dealii::VectorOperation::add);
  };

  set_value(*fixed_term_ptr, 1);
  set_value(*fission_term_ptr, 2);
  set_value(*other_term_ptr, 3);
  set_value(vector_1, 6);

  test_linear_term.SetFixedTermPtr({0,0}, fixed_term_ptr);
  test_linear_term.SetVariableTermPtr({0,0},
                                      VariableTerms::kFissionSource,
                                      fission_term_ptr);
  test_linear_term.SetVariableTermPtr({0,0},
                                      VariableTerms::kOther,
                                      other_term_ptr);

  auto term_vector_ptr = test_linear_term.GetFullTermPtr({0,0});
  EXPECT_TRUE(bart::testing::CompareMPIVectors(vector_1, *term_vector_ptr));
}

} // namespace