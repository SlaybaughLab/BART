#include "results/output_dealii_vtu.h"

namespace bart {

namespace results {

template <int dim>
OutputDealiiVtu<dim>::OutputDealiiVtu(
    const std::shared_ptr<domain::DefinitionI<dim>> &domain_ptr)
    : domain_ptr_(domain_ptr) {}

template <int dim>
void OutputDealiiVtu<dim>::AddData(bart::system::System &to_output) {
  //data_out_.attach_dof_handler(domain_ptr_->dof_handler());
}

template <int dim>
void OutputDealiiVtu<dim>::WriteData(std::ostream &output_stream) const {

}

template class OutputDealiiVtu<1>;
template class OutputDealiiVtu<2>;
template class OutputDealiiVtu<3>;

} // namespace results

} // namespace bart
