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
};

TYPED_TEST_SUITE(FormulationFactoryTests, bart::testing::AllDimensions);

TYPED_TEST(FormulationFactoryTests, MakeSAAFFormulationPtr) {
  constexpr int dim = this->dim;
  using FiniteElementType = NiceMock<domain::finite_element::FiniteElementMock<dim>>;
  using MaterialType = NiceMock<btest::MockMaterial>;
  using QuadratureSetType = NiceMock<quadrature::QuadratureSetMock<dim>>;

  auto finite_element_ptr = std::make_shared<FiniteElementType>();
  auto material_ptr = std::make_shared<MaterialType>();
  auto cross_section_ptr = std::make_shared<data::CrossSections>(*material_ptr);
  auto quadrature_ptr = std::make_shared<QuadratureSetType>();

  using BaseType = formulation::angular::SelfAdjointAngularFluxI<dim>;
  using ExpectedType = formulation::angular::SelfAdjointAngularFlux<dim>;

  std::unique_ptr<BaseType> returned_ptr = nullptr;
  EXPECT_NO_THROW({
    returned_ptr = std::move(formulation::factory::MakeSAAFFormulationPtr<dim>(
        finite_element_ptr,
        cross_section_ptr,
        quadrature_ptr));
  });
  ASSERT_NE(returned_ptr, nullptr);
  ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(returned_ptr.get()));
}


} // namespace
