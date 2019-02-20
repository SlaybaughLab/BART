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
  // Get material mapping from input file
  std::string material_mapping =
      GetMaterialMapFromFile(parameters_->MaterialMapFilename());
  // Instantiate mesh
  std::unique_ptr<domain::MeshI<dim>> mesh =
      std::make_unique<domain::MeshCartesian<dim>>(parameters_->SpatialMax(),
                                                   parameters_->NCells(),
                                                   material_mapping);
  
  // Instantiate FiniteElement
  std::unique_ptr<domain::FiniteElementI<dim>> finite_element =
      std::make_unique<domain::FiniteElementGaussian<dim>>(
          parameters_->Discretization(),
          parameters_->FEPolynomialDegree());
  
  // Instantiate domain definition
  auto definition = std::make_unique<domain::Definition<dim>>(
      mesh, finite_element);
  return definition;
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
