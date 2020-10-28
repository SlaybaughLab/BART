#ifndef BART_SRC_DOMAIN_MESH_FACTORY_HPP_
#define BART_SRC_DOMAIN_MESH_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.h"

namespace bart {

namespace domain {

namespace mesh {

template <int dim>
class MeshI;

enum class MeshName {
  kCartesian = 0,
};

template <int dim, typename ...T>
class MeshIFactory
    : public utility::factory::AutoRegisteringFactory<MeshName, std::unique_ptr<MeshI<dim>>(*)(T...)> {};

} // namespace mesh

} // namespace domain

} // namespace bart

#endif //BART_SRC_DOMAIN_MESH_FACTORY_HPP_
