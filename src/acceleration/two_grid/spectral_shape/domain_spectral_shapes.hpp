#ifndef BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_DOMAIN_SPECTRAL_SHAPES_HPP_
#define BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_DOMAIN_SPECTRAL_SHAPES_HPP_

#include "acceleration/two_grid/spectral_shape/domain_spectral_shapes_i.hpp"

namespace bart::acceleration::two_grid::spectral_shape {

/*! \brief Default implementation for domain spectral shapes calculator */
template <int dim>
class DomainSpectralShapes : public DomainSpectralShapesI<dim> {
 public:
  using typename DomainSpectralShapesI<dim>::Domain;
  using typename DomainSpectralShapesI<dim>::MaterialToGroupSpectralShapeMap;
  using typename DomainSpectralShapesI<dim>::GroupToDomainSpectralShapeMap;

  auto CalculateDomainSpectralShapes(const MaterialToGroupSpectralShapeMap&,
                                     const Domain&) const -> GroupToDomainSpectralShapeMap override;
};

} // namespace bart::acceleration::two_grid::spectral_shape

#endif //BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_DOMAIN_SPECTRAL_SHAPES_HPP_
