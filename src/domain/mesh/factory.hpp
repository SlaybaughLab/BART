#ifndef BART_SRC_DOMAIN_MESH_FACTORY_HPP_
#define BART_SRC_DOMAIN_MESH_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.hpp"

namespace bart::domain::mesh {

template <int dim>
class MeshI;

enum class MeshName {
  kCartesian = 0,
};

template <int dim, typename ...T>
class MeshIFactory
    : public utility::factory::AutoRegisteringFactory<MeshName, std::unique_ptr<MeshI<dim>>(*)(T...)> {};

[[nodiscard]] inline auto to_string(MeshName to_convert) -> std::string {
  switch (to_convert) {
    case (MeshName::kCartesian):
      return std::string{"MeshName::kCartesian"};
  }
  return std::string{"Unknown MeshName conversion to string requested."};
}

} // namespace bart::domain::mesh

#endif //BART_SRC_DOMAIN_MESH_FACTORY_HPP_
