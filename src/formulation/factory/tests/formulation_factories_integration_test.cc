#include "formulation/factory/formulation_factories.h"

#include "test_helpers/gmock_wrapper.h"

// Built by factory
#include "formulation/angular/self_adjoint_angular_flux.h"

// Dependencies and mocks
#include "data/cross_sections.h"
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
  using FiniteElementType = NiceMock<domain::finite_element::FiniteElementMock<dim>>;
  using MaterialType = NiceMock<btest::MockMaterial>;
  using QuadratureSetType = NiceMock<quadrature::QuadratureSetMock<dim>>;

  std::shared_ptr<FiniteElementType> finite_element_ptr_;
  std::shared_ptr<MaterialType> material_ptr_;
  std::shared_ptr<data::CrossSections> cross_section_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationFactoryTests<DimensionWrapper>::SetUp() {
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


} // namespace
