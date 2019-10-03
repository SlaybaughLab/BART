#ifndef BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_
#define BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_

namespace bart {

namespace quadrature {

template <int dim>
using CartesianPosition = utility::NamedType<std::array<double, dim>,
    struct CartesianPositionParameter>;

using Weight = utility::NamedType<double, struct WeightParameter>;

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_QUADRATURE_TYPES_H_
