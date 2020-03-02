#ifndef BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_
#define BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_

#include <memory>

#include "formulation/scalar/diffusion_i.h"
#include "formulation/stamper_i.h"
#include "formulation/updater/fixed_updater_i.h"
#include "formulation/updater/scattering_source_updater_i.h"
#include "formulation/updater/fission_source_updater_i.h"
#include "quadrature/quadrature_set_i.h"

namespace bart {

namespace formulation {

namespace updater {

template <int dim>
class DiffusionUpdater {
 public:
};

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_
