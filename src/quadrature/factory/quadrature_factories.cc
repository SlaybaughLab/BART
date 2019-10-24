#include "quadrature/factory/quadrature_factories.h"

#include "quadrature/ordinate.h"
#include "quadrature/quadrature_point.h"
#include "quadrature/quadrature_set.h"
#include "quadrature/angular/scalar_angular.h"
#include "quadrature/angular/level_symmetric_gaussian.h"
#include "quadrature/utility/quadrature_utilities.h"

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

template <int dim>
std::shared_ptr<QuadratureSetI<dim>> MakeQuadratureSetPtr(
    const QuadratureSetImpl type) {

  std::shared_ptr<QuadratureSetI<dim>> quadrature_set_ptr = nullptr;

  if (type == QuadratureSetImpl::kDefault) {
    quadrature_set_ptr = std::make_shared<QuadratureSet<dim>>();
  }

  return quadrature_set_ptr;
}

template <int dim>
void FillQuadratureSet(
    QuadratureSetI<dim>* to_fill,
    const std::vector<std::pair<quadrature::CartesianPosition<dim>, quadrature::Weight>>& point_vector) {

  for (const auto& [position, weight] : point_vector) {
    auto ordinate_ptr = MakeOrdinatePtr<dim>();
    ordinate_ptr->set_cartesian_position(position);

    auto quadrature_point_ptr = MakeQuadraturePointPtr<dim>();
    quadrature_point_ptr->SetOrdinate(ordinate_ptr).SetWeight(weight);

    to_fill->AddPoint(quadrature_point_ptr);

    auto reflection_ptr = MakeOrdinatePtr<dim>();
    reflection_ptr->set_cartesian_position(CartesianPosition<dim>(
        utility::ReflectAcrossOrigin<dim>(*ordinate_ptr)));
    auto reflection_point_ptr = MakeQuadraturePointPtr<dim>();
    reflection_point_ptr->SetOrdinate(reflection_ptr).SetWeight(weight);

    to_fill->AddPoint(reflection_point_ptr);

    to_fill->SetReflection(quadrature_point_ptr, reflection_point_ptr);
  }
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

template std::shared_ptr<QuadratureSetI<1>> MakeQuadratureSetPtr(const QuadratureSetImpl);
template std::shared_ptr<QuadratureSetI<2>> MakeQuadratureSetPtr(const QuadratureSetImpl);
template std::shared_ptr<QuadratureSetI<3>> MakeQuadratureSetPtr(const QuadratureSetImpl);

template void FillQuadratureSet<1>(QuadratureSetI<1>*, const std::vector<std::pair<quadrature::CartesianPosition<1>, quadrature::Weight>>&);
template void FillQuadratureSet<2>(QuadratureSetI<2>*, const std::vector<std::pair<quadrature::CartesianPosition<2>, quadrature::Weight>>&);
template void FillQuadratureSet<3>(QuadratureSetI<3>*, const std::vector<std::pair<quadrature::CartesianPosition<3>, quadrature::Weight>>&);

} // namespace factory

} // namespace quadrature

} // namespace bart