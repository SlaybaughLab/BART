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
  using typename FrameworkBuilderI<dim>::Stamper;

  MOCK_METHOD(std::unique_ptr<Domain>, BuildDomain, (const FrameworkParameters::DomainSize,
      const FrameworkParameters::NumberOfCells, const std::shared_ptr<FiniteElement>&,
      const std::string material_mapping), (override));
  MOCK_METHOD(std::unique_ptr<FiniteElement>, BuildFiniteElement, (const problem::CellFiniteElementType,
      const problem::DiscretizationType, const FrameworkParameters::PolynomialDegree), (override));
  MOCK_METHOD(std::shared_ptr<QuadratureSet>, BuildQuadratureSet, (const problem::AngularQuadType,
      const FrameworkParameters::AngularQuadratureOrder), (override));
  MOCK_METHOD(std::unique_ptr<Stamper>, BuildStamper, (const std::shared_ptr<Domain>&), (override));
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_BUILDER_MOCK_HPP_
