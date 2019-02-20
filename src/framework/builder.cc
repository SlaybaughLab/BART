#include "builder.h"

#include "../domain/mesh_cartesian.h"

namespace bart {

namespace framework {

template <int dim>
Builder<dim>::Builder(std::shared_ptr<problem::ParametersI> parameters)
    : parameters_(parameters) {}

template <int dim>
std::unique_ptr<domain::Definition<dim>> Builder<dim>::BuildDefinition() const {
  domain::MeshCartesian mesh(parameters_->SpatialMax(), parameters_->NCells());
  
  
}

} // namespace framework

} // namespace bart 
