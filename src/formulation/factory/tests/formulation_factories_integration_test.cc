#include "formulation/factory/formulation_factories.h"

#include "test_helpers/gmock_wrapper.h"

// Built by factory
#include "formulation/stamper.h"
#include "formulation/angular/self_adjoint_angular_flux.h"
#include "formulation/scalar/diffusion.h"


// Dependencies and mocks
#include "data/cross_sections.h"
#include "domain/tests/definition_mock.h"
#include "domain/finite_element/tests/finite_element_mock.h"
#include "material/tests/mock_material.h"
#include "quadrature/tests/quadrature_set_mock.h"


namespace  {

using namespace bart;

using ::testing::NiceMock;

template <typename DimensionWrapper>
class FormulationFactoryTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using DefinitionType = NiceMock<domain::DefinitionMock<dim>>;
  using FiniteElementType = NiceMock<domain::finite_element::FiniteElementMock<dim>>;
  using MaterialType = NiceMock<btest::MockMaterial>;
  using QuadratureSetType = NiceMock<quadrature::QuadratureSetMock<dim>>;

  std::shared_ptr<DefinitionType> definition_ptr_;
  std::shared_ptr<FiniteElementType> finite_element_ptr_;
  std::shared_ptr<MaterialType> material_ptr_;
  std::shared_ptr<data::CrossSections> cross_section_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationFactoryTests<DimensionWrapper>::SetUp() {
  definition_ptr_ = std::make_shared<DefinitionType>();
  finite_element_ptr_ = std::make_shared<FiniteElementType>();
  material_ptr_ = std::make_shared<MaterialType>();
  cross_section_ptr_ = std::make_shared<data::CrossSections>(*material_ptr_);
  quadrature_ptr_ = std::make_shared<QuadratureSetType>();
}

TYPED_TEST_SUITE(FormulationFactoryTests, bart::testing::AllDimensions);

TYPED_TEST(FormulationFactoryTests, MakeSAAFFormulationPtr) {
  constexpr int dim = this->dim;

  using BaseType = formulation::angular::SelfAdjointAngularFluxI<dim>;
  using ExpectedType = formulation::angular::SelfAdjointAngularFlux<dim>;

  std::unique_ptr<BaseType> returned_ptr = nullptr;
  EXPECT_NO_THROW({
    returned_ptr = std::move(formulation::factory::MakeSAAFFormulationPtr<dim>(
        this->finite_element_ptr_,
        this->cross_section_ptr_,
        this->quadrature_ptr_));
  });
  ASSERT_NE(returned_ptr, nullptr);
  ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(returned_ptr.get()));
}

TYPED_TEST(FormulationFactoryTests, MakeDiffusionPtr) {
  constexpr int dim = this->dim;
  using BaseType = formulation::scalar::DiffusionI<dim>;
  using ExpectedType = formulation::scalar::Diffusion<dim>;

  std::unique_ptr<BaseType> returned_ptr = nullptr;
  EXPECT_NO_THROW({
    returned_ptr = std::move(formulation::factory::MakeDiffusionPtr<dim>(
        this->finite_element_ptr_,
        this->cross_section_ptr_));
  });
  ASSERT_NE(returned_ptr, nullptr);
  ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(returned_ptr.get()));
}

TYPED_TEST(FormulationFactoryTests, MakeStamperPtr) {
  constexpr int dim = this->dim;
  using BaseType = formulation::StamperI<dim>;
  using ExpectedType = formulation::Stamper<dim>;

  std::unique_ptr<BaseType> returned_ptr = nullptr;
  EXPECT_NO_THROW({
    returned_ptr = std::move(formulation::factory::MakeStamperPtr<dim>(
        this->definition_ptr_));
  });
  ASSERT_NE(returned_ptr, nullptr);
  ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(returned_ptr.get()));
}

} // namespace
