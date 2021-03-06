#include "instrumentation/outstream/vector_to_vtu.hpp"

#include <filesystem>

namespace bart::instrumentation::outstream {

template<int dim>
VectorToVTU<dim>::VectorToVTU(std::shared_ptr<Definition> definition_ptr,
                              std::string data_name,
                              std::string directory,
                              std::string filename_base)
    : definition_ptr_(definition_ptr), data_name_(data_name), directory_(directory), filename_base_(filename_base) {}

template <int dim>
bool VectorToVTU<dim>::is_registered_ = Factory::get()
    .RegisterConstructor(OutstreamName::kVectorToVTU, [](std::shared_ptr<Definition> definition_ptr,
                                                         std::string data_name, std::string directory,
                                                         std::string filename_base)
                                                         -> std::unique_ptr<OutstreamI<Vector>> {
                           return std::make_unique<VectorToVTU<dim>>(definition_ptr,
                                                                     data_name, directory, filename_base);
    });

template <int dim>
auto VectorToVTU<dim>::Output(const Vector& to_output) -> VectorToVTU& {
  dealii::DataOut<dim> data_out;

  data_out.attach_dof_handler(definition_ptr_->dof_handler());
  data_out.add_data_vector(to_output, data_name_);
  data_out.build_patches();

  if (!std::filesystem::exists(directory_)) {
    std::filesystem::create_directories(directory_);
  }
  data_out.write_vtu_with_pvtu_record(directory_ + "/", filename_base_, counter_, MPI_COMM_WORLD, 4);
  ++counter_;
  return *this;
}

template class VectorToVTU<1>;
template class VectorToVTU<2>;
template class VectorToVTU<3>;

} // namespace bart::instrumentation::outstream

