#ifndef BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_
#define BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_

#include <array>

#include "utility/named_type.h"

namespace bart {

namespace quadrature {

template <int dim>
using CartesianPosition = bart::utility::NamedType<std::array<double, dim>,
    struct CartesianPositionParameter>;

using Weight = bart::utility::NamedType<double, struct WeightParameter>;

using Order = bart::utility::NamedType<int, struct OrderParameter>;

using PointIndex = bart::utility::NamedType<int, struct PointIndexParameter>;

using FillAllQuadrants =
    bart::utility::NamedType<bool, struct FillAllQuadrantsParameter>;

using QuadraturePointIndex =
    bart::utility::NamedType<int, struct QuadraturePointIndexParameter>;

enum class AngularQuadratureSetType {
  kNone = 0,
  kLevelSymmetricGaussian = 1,
  kGaussLegendre = 2,
};

enum class OrdinateType {
  kDefault = 0,
};

enum class QuadraturePointImpl {
  kDefault = 0,
};

enum class QuadratureSetImpl {
  kDefault = 0,
};

enum class MomentCalculatorImpl {
  kScalarMoment = 0,
  kZerothMomentOnly = 1,
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_
