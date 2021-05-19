#include "instrumentation/outstream/vector_map_to_vtu.hpp"

#include <filesystem>

namespace bart::instrumentation::outstream {

template<int dim>
VectorMapToVTU<dim>::VectorMapToVTU(std::shared_ptr<Definition> definition_ptr,
                                    const std::string data_name,
                                    const std::string directory,
                                    const std::string filename_base)
    : definition_ptr_(definition_ptr), data_name_(data_name), directory_(directory), filename_base_(filename_base) {}

template<int dim>
auto VectorMapToVTU<dim>::Output(const VectorMap& to_output) -> VectorMapToVTU& {
  dealii::DataOut<dim> data_out;

  data_out.attach_dof_handler(definition_ptr_->dof_handler());

  for (const auto& [group, vector] : to_output) {
    std::string vector_name{ data_name_ + "_group_" + std::to_string(group) };
    data_out.add_data_vector(vector, vector_name);
  }
  data_out.build_patches();

  if (!std::filesystem::exists(directory_)) {
    std::filesystem::create_directories(directory_);
  }
  data_out.write_vtu_with_pvtu_record(directory_ + "/", filename_base_, counter_, MPI_COMM_WORLD, 4);
  ++counter_;
  return *this;
}

template class VectorMapToVTU<1>;
template class VectorMapToVTU<2>;
template class VectorMapToVTU<3>;

} // namespace bart::instrumentation::outstream
