#include "bart_builder.h"

#include <unordered_map>

namespace bbuilders {

template <int dim>
void BuildMesh (dealii::ParameterHandler &prm,
    std::unique_ptr<MeshGenerator<dim>> &msh_ptr) {
  msh_ptr = std::unique_ptr<MeshGenerator<dim>> (new MeshGenerator<dim>(prm));
}

template <int dim>
std::unique_ptr<MeshGenerator<dim>> BuildMesh (
    dealii::ParameterHandler &prm) {
  return std::unique_ptr<MeshGenerator<dim>> (new MeshGenerator<dim>(prm));
}

}

// explicit instantiation for BuildMesh
template void bbuilders::BuildMesh<1> (dealii::ParameterHandler&,
    std::unique_ptr<MeshGenerator<1>> &);
template void bbuilders::BuildMesh<2> (dealii::ParameterHandler&,
    std::unique_ptr<MeshGenerator<2>> &);
template void bbuilders::BuildMesh<3> (dealii::ParameterHandler&,
    std::unique_ptr<MeshGenerator<3>> &);

template std::unique_ptr<MeshGenerator<1>> bbuilders::BuildMesh<1> (
    dealii::ParameterHandler&);
template std::unique_ptr<MeshGenerator<2>> bbuilders::BuildMesh<2> (
    dealii::ParameterHandler&);
template std::unique_ptr<MeshGenerator<3>> bbuilders::BuildMesh<3> (
    dealii::ParameterHandler&);
