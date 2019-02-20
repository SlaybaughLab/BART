#include "builder.h"

#include <fstream>

#include "../domain/mesh_cartesian.h"

namespace bart {

namespace framework {

template <int dim>
Builder<dim>::Builder(std::shared_ptr<problem::ParametersI> parameters)
    : parameters_(parameters) {}

template <int dim>
std::unique_ptr<domain::Definition<dim>> Builder<dim>::BuildDefinition() const {
  // Get material mapping from input file and build mesh
  std::string material_mapping =
      GetMaterialMapFromFile(parameters_->MaterialMapFilename());
  domain::MeshCartesian mesh(parameters_->SpatialMax(),
                             parameters_->NCells(),
                             material_mapping);

  
  
}

template <int dim>
std::string Builder<dim>::GetMaterialMapFromFile(std::string filename) const {
  std::ifstream input(filename);
  std::string return_string(std::istreambuf_iterator<char>(input),
                            std::istreambuf_iterator<char>());
  return return_string;
}

} // namespace framework

} // namespace bart 
