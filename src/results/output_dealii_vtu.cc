#include "results/output_dealii_vtu.h"

namespace bart {

namespace results {

template <int dim>
void OutputDealiiVtu<dim>::AddData(bart::system::System &to_output) {

}

template <int dim>
void OutputDealiiVtu<dim>::WriteData(std::ostream &output_stream) const {

}

template class OutputDealiiVtu<1>;
template class OutputDealiiVtu<2>;
template class OutputDealiiVtu<3>;

} // namespace results

} // namespace bart
