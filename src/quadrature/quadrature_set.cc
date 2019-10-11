#include "quadrature/quadrature_set.h"

namespace bart {

namespace quadrature {

template<int dim>
bool QuadratureSet<dim>::AddPoint(
    std::shared_ptr<QuadraturePointI<dim>> new_point_ptr) {
  auto status = quadrature_point_ptrs_.insert(new_point_ptr);
  return status.second;
}

template class QuadratureSet<1>;
template class QuadratureSet<2>;
template class QuadratureSet<3>;

} // namespace quadrature

} // namespace bart

