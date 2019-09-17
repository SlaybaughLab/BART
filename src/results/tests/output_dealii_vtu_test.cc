#include "results/output_dealii_vtu.h"

#include <memory>

#include "domain/tests/definition_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

template <typename DimensionWrapper>
class ResultsOutputDealiiVtuTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  std::shared_ptr<domain::DefinitionI<dim>> domain_ptr_;
  std::unique_ptr<results::OutputDealiiVtu<dim>> test_output_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void ResultsOutputDealiiVtuTest<DimensionWrapper>::SetUp() {
  domain_ptr_ = std::make_shared<domain::DefinitionMock<dim>>();

  test_output_ = std::make_unique<results::OutputDealiiVtu<dim>>(domain_ptr_);
}

TYPED_TEST_CASE(ResultsOutputDealiiVtuTest, bart::testing::AllDimensions);

TYPED_TEST(ResultsOutputDealiiVtuTest, Constructor) {
  constexpr int dim = this->dim;
  EXPECT_EQ(this->domain_ptr_.use_count(), 2);

  auto domain_ptr = this->test_output_->domain_ptr();

  ASSERT_NE(domain_ptr, nullptr);
  EXPECT_NE(nullptr, dynamic_cast<domain::DefinitionMock<dim>*>(domain_ptr));
}

}