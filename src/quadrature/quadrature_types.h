#ifndef BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_
#define BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_

#include <array>

#include "utility/named_type.h"

namespace bart {

namespace quadrature {

template <int dim>
using CartesianPosition = utility::NamedType<std::array<double, dim>,
    struct CartesianPositionParameter>;

using Weight = utility::NamedType<double, struct WeightParameter>;

using Order = utility::NamedType<int, struct OrderParameter>;

using PointIndex = utility::NamedType<int, struct PointIndexParameter>;

enum class AngularQuadratureSetType {
  kNone = 0,
  kScalar = 1,
  kLevelSymmetricGaussian = 2,
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

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_
