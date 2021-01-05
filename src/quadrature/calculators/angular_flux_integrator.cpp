#include "quadrature/calculators/angular_flux_integrator.hpp"

namespace bart::quadrature::calculators {

template<int dim>
AngularFluxIntegrator<dim>::AngularFluxIntegrator(std::shared_ptr<QuadratureSet> quadrature_set_ptr)
    : quadrature_set_ptr_(quadrature_set_ptr) {
  this->AssertPointerNotNull(quadrature_set_ptr_.get(), "quadrature set", "DriftDiffusionIntegratedFlux constructor");
}

template<int dim>
auto AngularFluxIntegrator<dim>::Integrate(const VectorMap& vector_map) const -> Vector {
  using Index = quadrature::QuadraturePointIndex;
  const VectorMap::size_type n_quadrature_points{ this->quadrature_set_ptr_->size() };
  AssertThrow(vector_map.size() == n_quadrature_points, dealii::ExcMessage("Error in DriftDiffusionIntegratedFlux "
                                                                           "integration, angular flux map is not the "
                                                                           "same size as the quadrature set"));
  const Vector::size_type vector_size{ vector_map.begin()->second->size()};
  Vector result_vector(vector_size);

  for (VectorMap::size_type i = 0; i < n_quadrature_points; ++i) {
    auto& quadrature_point = *quadrature_set_ptr_->GetQuadraturePoint(Index(i));
    const double weight{ quadrature_point.weight() };
    const auto position{ quadrature_point.cartesian_position_tensor() };
    result_vector.add(weight * position * position, *vector_map.at(Index(i)));
  }

  return result_vector;
}

template class AngularFluxIntegrator<1>;
template class AngularFluxIntegrator<2>;
template class AngularFluxIntegrator<3>;

} // namespace bart::quadrature::calculators
