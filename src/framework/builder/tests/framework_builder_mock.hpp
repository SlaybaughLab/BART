#ifndef BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_BUILDER_MOCK_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_BUILDER_MOCK_HPP_

#include "framework/builder/framework_builder_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::framework::builder {

template <int dim>
class FrameworkBuilderMock : public FrameworkBuilderI<dim> {
 public:
  using typename FrameworkBuilderI<dim>::Domain;
  using typename FrameworkBuilderI<dim>::FiniteElement;
  using typename FrameworkBuilderI<dim>::FrameworkI;
  using typename FrameworkBuilderI<dim>::QuadratureSet;
  using typename FrameworkBuilderI<dim>::SAAFFormulation;
  using typename FrameworkBuilderI<dim>::Stamper;

  using typename FrameworkBuilderI<dim>::AngularFluxStorage;

  using typename FrameworkBuilderI<dim>::UpdaterPointers;
  using typename FrameworkBuilderI<dim>::BoundaryConditionsUpdater;
  using typename FrameworkBuilderI<dim>::FissionSourceUpdater;
  using typename FrameworkBuilderI<dim>::FixedTermUpdater ;
  using typename FrameworkBuilderI<dim>::ScatteringSourceUpdater;

  MOCK_METHOD(std::unique_ptr<Domain>, BuildDomain, (const FrameworkParameters::DomainSize,
      const FrameworkParameters::NumberOfCells, const std::shared_ptr<FiniteElement>&,
      const std::string material_mapping), (override));
  MOCK_METHOD(std::unique_ptr<FiniteElement>, BuildFiniteElement, (const problem::CellFiniteElementType,
      const problem::DiscretizationType, const FrameworkParameters::PolynomialDegree), (override));
  MOCK_METHOD(std::shared_ptr<QuadratureSet>, BuildQuadratureSet, (const problem::AngularQuadType,
      const FrameworkParameters::AngularQuadratureOrder), (override));
  MOCK_METHOD(std::unique_ptr<SAAFFormulation>, BuildSAAFFormulation, (const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::CrossSections>&, const std::shared_ptr<QuadratureSet>&,
      const formulation::SAAFFormulationImpl), (override));
  MOCK_METHOD(std::unique_ptr<Stamper>, BuildStamper, (const std::shared_ptr<Domain>&), (override));
  MOCK_METHOD(UpdaterPointers, BuildUpdaterPointers, (std::unique_ptr<SAAFFormulation>,
      std::unique_ptr<Stamper>, const std::shared_ptr<QuadratureSet>&), (override));
  MOCK_METHOD(UpdaterPointers, BuildUpdaterPointers, (std::unique_ptr<SAAFFormulation>, std::unique_ptr<Stamper>,
      const std::shared_ptr<QuadratureSet>&, (const std::map<problem::Boundary, bool>)& reflective_boundaries,
      const AngularFluxStorage&), (override));
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_BUILDER_MOCK_HPP_
