
#include "quadrature/factory/quadrature_factories.h"

#include "quadrature/ordinate.h"
#include "quadrature/quadrature_point.h"
#include "quadrature/quadrature_set.h"
#include "quadrature/angular/level_symmetric_gaussian.h"
#include "quadrature/calculators/scalar_moment.h"
#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"
#include "quadrature/utility/quadrature_utilities.h"

namespace bart {

namespace quadrature {

namespace factory {

std::string unsupported_quadrature_error = "Error in MakeAngularQuadratureGeneratorPtr, unsupported type of AngularQuadrature<dim> requested";

template <int dim>
std::shared_ptr<OrdinateI<dim>> MakeOrdinatePtr(const OrdinateType type) {
  std::shared_ptr<OrdinateI<dim>> ordinate_ptr = nullptr;

  if (type == OrdinateType::kDefault) {
    ordinate_ptr = std::make_shared<Ordinate<dim>>();
  }

  return ordinate_ptr;
}

template <int dim>
std::shared_ptr<QuadraturePointI<dim>> MakeQuadraturePointPtr(
    const QuadraturePointImpl impl) {
  std::shared_ptr<QuadraturePointI<dim>> quadrature_point_ptr = nullptr;

  if (impl == QuadraturePointImpl::kDefault) {
    quadrature_point_ptr = std::make_shared<QuadraturePoint<dim>>();
  }

  return quadrature_point_ptr;
}

template <>
std::shared_ptr<QuadratureGeneratorI<3>> MakeAngularQuadratureGeneratorPtr(
    const Order order,
    const AngularQuadratureSetType type) {
  std::shared_ptr<QuadratureGeneratorI<3>> generator_ptr = nullptr;

  if (type == AngularQuadratureSetType::kLevelSymmetricGaussian) {
    generator_ptr = std::make_shared<angular::LevelSymmetricGaussian>(order);
  } else {
    AssertThrow(false, dealii::ExcMessage(unsupported_quadrature_error));
  }

  return generator_ptr;
}

template <int dim>
std::shared_ptr<QuadratureGeneratorI<dim>> MakeAngularQuadratureGeneratorPtr(
    const Order,
    const AngularQuadratureSetType) {
  std::shared_ptr<QuadratureGeneratorI<dim>> generator_ptr = nullptr;

  AssertThrow(false, dealii::ExcMessage(unsupported_quadrature_error));

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

template<int dim>
std::unique_ptr<calculators::SphericalHarmonicMomentsI> MakeMomentCalculator(
    const MomentCalculatorImpl impl,
    std::shared_ptr<QuadratureSetI<dim>> quadrature_set_ptr) {

  std::unique_ptr<calculators::SphericalHarmonicMomentsI> return_pointer =
      nullptr;

  if (impl == MomentCalculatorImpl::kScalarMoment) {
    return_pointer = std::move(
        std::make_unique<calculators::ScalarMoment>());
  } else if (impl == MomentCalculatorImpl::kZerothMomentOnly) {
    AssertThrow(quadrature_set_ptr != nullptr,
                dealii::ExcMessage("Error in factory building moment calculator, "
                                   "implementation requires quadrature set but provided"
                                   " set pointer is a nullptr"))
    return_pointer = std::move(
        std::make_unique<calculators::SphericalHarmonicZerothMoment<dim>>(
            quadrature_set_ptr));
  }

  return std::move(return_pointer);
}

template <int dim>
void FillQuadratureSet(
    QuadratureSetI<dim>* to_fill,
    const std::vector<std::pair<quadrature::CartesianPosition<dim>, quadrature::Weight>>& point_vector) {

  AssertThrow(to_fill->size() == 0,
      dealii::ExcMessage("Error in FillQuadratureSet, set is not empty"));

  for (const auto& [position, weight] : point_vector) {
    auto ordinate_ptr = MakeOrdinatePtr<dim>();
    ordinate_ptr->set_cartesian_position(position);

    auto quadrature_point_ptr = MakeQuadraturePointPtr<dim>();
    quadrature_point_ptr->SetOrdinate(ordinate_ptr).SetWeight(weight);

    to_fill->AddPoint(quadrature_point_ptr);

    auto reflection_ptr = MakeOrdinatePtr<dim>();
    reflection_ptr->set_cartesian_position(CartesianPosition<dim>(
        utility::ReflectAcrossOrigin<dim>(*ordinate_ptr)));

    // Only add if the point isn't its own reflection
    if (reflection_ptr->cartesian_position() !=
        ordinate_ptr->cartesian_position()) {
      auto reflection_point_ptr = MakeQuadraturePointPtr<dim>();
      reflection_point_ptr->SetOrdinate(reflection_ptr).SetWeight(weight);

      to_fill->AddPoint(reflection_point_ptr);

      to_fill->SetReflection(quadrature_point_ptr, reflection_point_ptr);
    }
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

template std::unique_ptr<calculators::SphericalHarmonicMomentsI> MakeMomentCalculator<1>(const MomentCalculatorImpl, std::shared_ptr<QuadratureSetI<1>>);
template std::unique_ptr<calculators::SphericalHarmonicMomentsI> MakeMomentCalculator<2>(const MomentCalculatorImpl, std::shared_ptr<QuadratureSetI<2>>);
template std::unique_ptr<calculators::SphericalHarmonicMomentsI> MakeMomentCalculator<3>(const MomentCalculatorImpl, std::shared_ptr<QuadratureSetI<3>>);

template void FillQuadratureSet<1>(QuadratureSetI<1>*, const std::vector<std::pair<quadrature::CartesianPosition<1>, quadrature::Weight>>&);
template void FillQuadratureSet<2>(QuadratureSetI<2>*, const std::vector<std::pair<quadrature::CartesianPosition<2>, quadrature::Weight>>&);
template void FillQuadratureSet<3>(QuadratureSetI<3>*, const std::vector<std::pair<quadrature::CartesianPosition<3>, quadrature::Weight>>&);

} // namespace factory

} // namespace quadrature

} // namespace bart