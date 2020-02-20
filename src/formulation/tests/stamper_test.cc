#include "formulation/stamper.h"

#include "domain/tests/definition_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

template <typename DimensionWrapper>
class FormulationStamperTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using DomainDefinitionType = domain::DefinitionMock<dim>;

  std::shared_ptr<DomainDefinitionType> domain_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationStamperTest<DimensionWrapper>::SetUp() {
  domain_ptr_ = std::make_shared<DomainDefinitionType>();
}

TYPED_TEST_SUITE(FormulationStamperTest, bart::testing::AllDimensions);

TYPED_TEST(FormulationStamperTest, Constructor) {
  using StamperType = formulation::Stamper<this->dim>;
  std::shared_ptr<StamperType> stamper_ptr;
  EXPECT_NO_THROW({
    stamper_ptr = std::make_shared<StamperType>(this->domain_ptr_);
  });
  ASSERT_NE(stamper_ptr->domain_ptr(), nullptr);
  EXPECT_EQ(stamper_ptr->domain_ptr(), this->domain_ptr_.get());
}

TYPED_TEST(FormulationStamperTest, BadDependency) {
  using StamperType = formulation::Stamper<this->dim>;
  std::shared_ptr<StamperType> stamper_ptr;
  EXPECT_ANY_THROW({
    stamper_ptr = std::make_shared<StamperType>(nullptr);
  });
}

} // namespace
