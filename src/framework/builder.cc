#include "builder.h"

#include <fstream>

#include "../domain/finite_element_gaussian.h"
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
  auto mesh = std::make_unique<domain::MeshCartesian<dim>>(
      parameters_->SpatialMax(), parameters_->NCells(), material_mapping);
  auto finite_element = std::make_unique<domain::FiniteElementGaussian<dim>>(
      parameters_->Discretization(), parameters_->FEPolynomialDegree());

}

template <int dim>
std::string Builder<dim>::GetMaterialMapFromFile(std::string filename) const {
  std::ifstream input(filename);
  std::string return_string((std::istreambuf_iterator<char>(input)),
                            (std::istreambuf_iterator<char>()));
  return return_string;
}

template class Builder<1>;
template class Builder<2>;

} // namespace framework

} // namespace bart 
