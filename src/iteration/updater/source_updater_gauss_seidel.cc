#include "iteration/updater/source_updater_gauss_seidel.h"

namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
void SourceUpdaterGaussSeidel<StamperType>::UpdateScatteringSource(System& system) {

}

template class SourceUpdaterGaussSeidel<formulation::CFEMStamperI>;

} // namespace updater

} // namespace iteration

} // namespace bart