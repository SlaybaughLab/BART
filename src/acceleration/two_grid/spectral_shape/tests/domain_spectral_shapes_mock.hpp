#ifndef BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_TESTS_DOMAIN_SPECTRAL_SHAPES_MOCK_HPP_
#define BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_TESTS_DOMAIN_SPECTRAL_SHAPES_MOCK_HPP_

#include "acceleration/two_grid/spectral_shape/domain_spectral_shapes_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::acceleration::two_grid::spectral_shape {

template <int dim>
class DomainSpectralShapesMock : public DomainSpectralShapesI<dim> {
 public:
  using typename DomainSpectralShapesI<dim>::Domain;
  using typename DomainSpectralShapesI<dim>::MaterialToGroupSpectralShapeMap;
  using typename DomainSpectralShapesI<dim>::GroupToDomainSpectralShapeMap;

  MOCK_METHOD(GroupToDomainSpectralShapeMap, CalculateDomainSpectralShapes, (const MaterialToGroupSpectralShapeMap&,
      const Domain&), (const, override));
};

} // namespace bart::acceleration::two_grid::spectral_shape

#endif //BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_TESTS_DOMAIN_SPECTRAL_SHAPES_MOCK_HPP_
