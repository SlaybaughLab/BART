#include <ostream>

#include "system/system.hpp"
#include "results/output_dealii_vtu.h"

namespace bart {

namespace results {

template <int dim>
OutputDealiiVtu<dim>::OutputDealiiVtu(
    const std::shared_ptr<domain::DefinitionI<dim>> &domain_ptr)
    : domain_ptr_(domain_ptr) {}

template <int dim>
void OutputDealiiVtu<dim>::AddData(bart::system::System &to_output) {
  data_out_.attach_dof_handler(domain_ptr_->dof_handler());

  for (int group = 0; group < to_output.total_groups; ++group) {
    std::ostringstream group_flux_label;
    group_flux_label << "scalar_flux_group_" << group;
    system::moments::MomentIndex index{group, 0, 0};
    auto& moment = to_output.current_moments->GetMoment(index);
    data_out_.add_data_vector(moment, group_flux_label.str());
  }

  data_out_.build_patches();
}

template <int dim>
void OutputDealiiVtu<dim>::WriteData(std::ostream &output_stream) const {
  data_out_.write_vtu(output_stream);
}

template<int dim>
void OutputDealiiVtu<dim>::WriteMasterFile(
    std::ostream &output_stream,
    std::vector<std::string> filenames) const {
  data_out_.write_pvtu_record(output_stream, filenames);
}

template class OutputDealiiVtu<1>;
template class OutputDealiiVtu<2>;
template class OutputDealiiVtu<3>;

} // namespace results

} // namespace bart
