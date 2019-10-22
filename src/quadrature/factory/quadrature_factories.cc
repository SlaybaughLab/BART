#include "quadrature/factory/quadrature_factories.h"

#include "quadrature/ordinate.h"
#include "quadrature/quadrature_point.h"
#include "quadrature/angular/scalar_angular.h"
#include "quadrature/angular/level_symmetric_gaussian.h"

namespace bart {

namespace quadrature {

namespace factory {

template <int dim>
std::shared_ptr<OrdinateI<dim>> MakeOrdinatePtr(const OrdinateType type) {

  std::shared_ptr<OrdinateI<dim>> ordinate_ptr = nullptr;

  if (type == OrdinateType::kDefault) {
    // Default implementation
    ordinate_ptr = std::make_shared<Ordinate<dim>>();
  }

  return ordinate_ptr;
}

template <int dim>
std::shared_ptr<QuadraturePointI<dim>> MakeQuadraturePointPtr(
    const QuadraturePointImpl impl) {
  std::shared_ptr<QuadraturePointI<dim>> quadrature_point_ptr = nullptr;

  if (impl == QuadraturePointImpl::kDefault) {
    // Default implementation
    quadrature_point_ptr = std::make_shared<QuadraturePoint<dim>>();
  }

  return quadrature_point_ptr;
}

template <>
std::shared_ptr<QuadratureGeneratorI<3>> MakeAngularQuadratureGeneratorPtr(
    const Order order,
    const AngularQuadratureSetType type) {
  std::shared_ptr<QuadratureGeneratorI<3>> generator_ptr = nullptr;

  if (type == AngularQuadratureSetType::kScalar) {
    generator_ptr = std::make_shared<angular::ScalarAngular<3>>();
  } else if (type == AngularQuadratureSetType::kLevelSymmetricGaussian) {
    generator_ptr = std::make_shared<angular::LevelSymmetricGaussian>(order);
  } else {
    AssertThrow(false,
                dealii::ExcMessage("Error in MakeAngularQuadratureGeneratorPtr, "
                                   "unsupported type of AngularQuadrature requested"));
  }

  return generator_ptr;
}

template <int dim>
std::shared_ptr<QuadratureGeneratorI<dim>> MakeAngularQuadratureGeneratorPtr(
    const Order order,
    const AngularQuadratureSetType type) {
  std::shared_ptr<QuadratureGeneratorI<dim>> generator_ptr = nullptr;

  if (type == AngularQuadratureSetType::kScalar) {
    generator_ptr = std::make_shared<angular::ScalarAngular<dim>>();
  } else {
    AssertThrow(false,
        dealii::ExcMessage("Error in MakeAngularQuadratureGeneratorPtr, "
                           "unsupported type of AngularQuadrature requested"));
  }

  return generator_ptr;
}


template std::shared_ptr<OrdinateI<1>> MakeOrdinatePtr(const OrdinateType);
template std::shared_ptr<OrdinateI<2>> MakeOrdinatePtr(const OrdinateType);
template std::shared_ptr<OrdinateI<3>> MakeOrdinatePtr(const OrdinateType);

template std::shared_ptr<QuadraturePointI<1>> MakeQuadraturePointPtr(const QuadraturePointImpl);
template std::shared_ptr<QuadraturePointI<2>> MakeQuadraturePointPtr(const QuadraturePointImpl);
template std::shared_ptr<QuadraturePointI<3>> MakeQuadraturePointPtr(const QuadraturePointImpl);

template std::shared_ptr<QuadratureGeneratorI<1>> MakeAngularQuadratureGeneratorPtr(const Order, const AngularQuadratureSetType);
template std::shared_ptr<QuadratureGeneratorI<2>> MakeAngularQuadratureGeneratorPtr(const Order, const AngularQuadratureSetType);
template std::shared_ptr<QuadratureGeneratorI<3>> MakeAngularQuadratureGeneratorPtr(const Order, const AngularQuadratureSetType);

} // namespace factory

} // namespace quadrature

} // namespace bart