#ifndef BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_DOMAIN_SPECTRAL_SHAPES_I_HPP_
#define BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_DOMAIN_SPECTRAL_SHAPES_I_HPP_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/vector.h>

#include "domain/domain_i.hpp"
#include "utility/has_description.h"

namespace bart::acceleration::two_grid::spectral_shape {

/*! \brief Interface for classes that calculate the group spectral shapes for a domain. */
template <int dim>
class DomainSpectralShapesI : public utility::HasDescription {
 public:
  using Domain = domain::DomainI<dim>;
  //! Mapping of material IDs to a vector of spectral shape values by group
  using MaterialToGroupSpectralShapeMap = std::unordered_map<int, std::vector<double>>;
  //! Mapping of groups to a vector of spectral shape values by global DOF
  using GroupToDomainSpectralShapeMap = std::unordered_map<int, dealii::Vector<double>>;
  virtual ~DomainSpectralShapesI() = default;
  virtual auto CalculateDomainSpectralShapes(const MaterialToGroupSpectralShapeMap&,
                                             const Domain&) const -> GroupToDomainSpectralShapeMap = 0;
};

} // namespace bart::acceleration::two_grid::spectral_shape

#endif //BART_SRC_ACCELERATION_TWO_GRID_SPECTRAL_SHAPE_DOMAIN_SPECTRAL_SHAPES_I_HPP_
