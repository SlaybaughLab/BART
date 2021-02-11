#include "formulation/factory/formulation_factories.h"

#include "test_helpers/gmock_wrapper.h"

// Built by factory
#include "formulation/stamper.h"
#include "formulation/angular/self_adjoint_angular_flux.h"
#include "formulation/scalar/diffusion.h"
#include "formulation/updater/saaf_updater.h"

// Dependencies and mocks
#include "data/cross_sections.h"
#include "domain/tests/definition_mock.h"
#include "domain/finite_element/tests/finite_element_mock.hpp"
#include "formulation/angular/tests/self_adjoint_angular_flux_mock.h"
#include "formulation/scalar/tests/diffusion_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "material/tests/material_mock.hpp"
#include "quadrature/tests/quadrature_set_mock.h"

namespace  {

using namespace bart;

using ::testing::NiceMock, ::testing::Return;

template <typename DimensionWrapper>
class FormulationFactoryTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using DefinitionType = NiceMock<domain::DefinitionMock<dim>>;
  using FiniteElementType = NiceMock<domain::finite_element::FiniteElementMock<dim>>;
  using DiffusionFormulationType = NiceMock<formulation::scalar::DiffusionMock<dim>>;
  using SAAFFormulationType = NiceMock<formulation::angular::SelfAdjointAngularFluxMock<dim>>;
  using MaterialType = NiceMock<btest::MockMaterial>;
  using QuadratureSetType = NiceMock<quadrature::QuadratureSetMock<dim>>;
  using StamperType = NiceMock<formulation::StamperMock<dim>>;

  std::shared_ptr<DefinitionType> definition_ptr_;
  std::shared_ptr<FiniteElementType> finite_element_ptr_;
  std::unique_ptr<DiffusionFormulationType> diffusion_formulation_ptr_;
  std::unique_ptr<SAAFFormulationType> saaf_formulation_ptr_;
  std::shared_ptr<MaterialType> material_ptr_;
  std::shared_ptr<data::CrossSections> cross_section_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_ptr_;
  std::unique_ptr<StamperType> stamper_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationFactoryTests<DimensionWrapper>::SetUp() {
  definition_ptr_ = std::make_shared<DefinitionType>();
  finite_element_ptr_ = std::make_shared<FiniteElementType>();
  diffusion_formulation_ptr_ = std::move(
      std::make_unique<DiffusionFormulationType>());
  saaf_formulation_ptr_ = std::move(std::make_unique<SAAFFormulationType>());
  material_ptr_ = std::make_shared<MaterialType>();
  cross_section_ptr_ = std::make_shared<data::CrossSections>(*material_ptr_);
  quadrature_ptr_ = std::make_shared<QuadratureSetType>();
  stamper_ptr_ = std::move(std::make_unique<StamperType>());
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

TYPED_TEST(FormulationFactoryTests, MakeDiffusionUpdaterPtr) {
  constexpr int dim = this->dim;
  using ExpectedType = formulation::updater::DiffusionUpdater<dim>;

  std::unique_ptr<ExpectedType> returned_ptr = nullptr;

  EXPECT_NO_THROW({
    returned_ptr = std::move(formulation::factory::MakeDiffusionUpdater<dim>(
        std::move(this->diffusion_formulation_ptr_),
        std::move(this->stamper_ptr_)));
  });
  ASSERT_NE(returned_ptr, nullptr);
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

TYPED_TEST(FormulationFactoryTests, MakeSAAFUpdater) {
  constexpr int dim = this->dim;
  using ExpectedType = formulation::updater::SAAFUpdater<dim>;

  std::unique_ptr<ExpectedType> returned_ptr = nullptr;
  EXPECT_NO_THROW({
    returned_ptr = std::move(formulation::factory::MakeSAAFUpdater<dim>(
        std::move(this->saaf_formulation_ptr_),
        std::move(this->stamper_ptr_),
        this->quadrature_ptr_));
  });
  ASSERT_NE(returned_ptr, nullptr);
}

} // namespace
