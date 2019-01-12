#include "bart_builder.h"

#include <unordered_map>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>

namespace bbuilders {
template <int dim>
void BuildFESpaces (const dealii::ParameterHandler &prm,
    std::unordered_map<std::string, dealii::FiniteElement<dim,dim>*>&fe_ptrs) {
  // getting parameter values
  const bool do_nda = prm.get_bool ("do nda");
  const int p_order = prm.get_integer("finite element polynomial degree");
  const std::string ho_discretization = prm.get ("ho spatial discretization");
  const std::string nda_discretization = prm.get ("nda spatial discretization");
  const std::string ho_equ_name = prm.get ("transport model");

  std::unordered_map<std::string, int> discretization_ind = {
      {"cfem", 0}, {"dfem", 1}, {"cmfd", 2}, {"rtk", 3}};

  switch (discretization_ind[ho_discretization]) {
    case 0:
      fe_ptrs[ho_equ_name] = new dealii::FE_Q<dim> (p_order);
      break;

    case 1:
      fe_ptrs[ho_equ_name] = new dealii::FE_DGQ<dim> (p_order);
      break;

    default:
      AssertThrow (false,
          dealii::ExcMessage("Invalid HO discretization name"));
      break;
  }

  if (do_nda) {
    switch (discretization_ind[nda_discretization]) {
      case 0:
        fe_ptrs["nda"] = new dealii::FE_Q<dim> (p_order);
        break;

      case 1:
        fe_ptrs["nda"] = new dealii::FE_DGQ<dim> (p_order);
        break;

      case 2:
        fe_ptrs["nda"] = new dealii::FE_DGQ<dim> (0);
        break;

      case 3:
        fe_ptrs["nda"] = new dealii::FE_RaviartThomas<dim> (p_order);
        break;

      default:
        AssertThrow (false,
            dealii::ExcMessage("Invalid NDA discretization name"));
        break;
    }
  }
}

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

// explicitly instantiate all builders using templates
// explicit instantiation for BuildFESpaces
template void bbuilders::BuildFESpaces<1> (const dealii::ParameterHandler&,
    std::unordered_map<std::string, dealii::FiniteElement<1, 1>*>&);
template void bbuilders::BuildFESpaces<2> (const dealii::ParameterHandler&,
    std::unordered_map<std::string, dealii::FiniteElement<2, 2>*>&);
template void bbuilders::BuildFESpaces<3> (const dealii::ParameterHandler&,
    std::unordered_map<std::string, dealii::FiniteElement<3, 3>*>&);

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
