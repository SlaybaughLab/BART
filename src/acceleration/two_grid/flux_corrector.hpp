#ifndef BART_SRC_ACCELERATION_FLUX_CORRECTOR_HPP_
#define BART_SRC_ACCELERATION_FLUX_CORRECTOR_HPP_

#include <unordered_map>

#include "acceleration/two_grid/flux_corrector_i.hpp"

//! Calculators for the two-grid acceleration method <a href="https://doi.org/10.13182/NSE115-253">Adams and Morel (2017)</a>
namespace bart::acceleration::two_grid {

class FluxCorrector : public FluxCorrectorI {
 public:
  using GroupToDomainSpectralShapeMap = std::unordered_map<int, dealii::Vector<double>>;
  explicit FluxCorrector(const GroupToDomainSpectralShapeMap& group_to_domain_spectral_shape_map)
      : group_to_domain_spectral_shape_map_(group_to_domain_spectral_shape_map) {}

  auto CorrectFlux(Vector &flux_to_correct, const Vector &error, const int group) const -> void override {
    auto scaled_error = error;
    scaled_error.scale(group_to_domain_spectral_shape_map_.at(group));
    flux_to_correct.add(1, scaled_error);
  }

  auto group_to_domain_spectral_shape_map() const { return group_to_domain_spectral_shape_map_; }
 private:
  GroupToDomainSpectralShapeMap group_to_domain_spectral_shape_map_;
};

} // namespace bart::acceleration::two_grid

#endif //BART_SRC_ACCELERATION_FLUX_CORRECTOR_HPP_
