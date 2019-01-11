#include "bart_builder.h"
#include "../aqdata/lsgc.h"
#include "../equation/even_parity.h"
#include "../equation/self_adjoint_angular_flux.h"
#include "../iteration/power_iteration.h"
#include "../iteration/gauss_seidel.h"
#include "../iteration/source_iteration.h"

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

template <int dim>
std::unordered_map<std::string, std::unique_ptr<EquationBase<dim>>> BuildEqu (
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr) {
  std::unordered_map<std::string, std::unique_ptr<EquationBase<dim>>> equ_ptrs;
  const std::string ho_equ_name(prm.get("transport model"));
  std::unordered_map<std::string, int> mp = {{"ep",0}, {"saaf",1}};
  switch (mp[ho_equ_name]) {
    case 0:
      equ_ptrs[ho_equ_name] = std::unique_ptr<EquationBase<dim>>(
          new EvenParity<dim>(ho_equ_name, prm, dat_ptr));
      break;
    case 1:
      equ_ptrs[ho_equ_name] = std::unique_ptr<EquationBase<dim>>(
          new SelfAdjointAngularFlux<dim>(ho_equ_name, prm, dat_ptr));
      break;
    default:
      break;
  }
  bool do_nda = prm.get_bool("do nda");
  if (do_nda) {
    AssertThrow(ho_equ_name!="diffusion",
        dealii::ExcMessage("NDA cannot be used with diffusion"));
    //TODO:
    // 1. Fill in NDA part once corresponding class is ready
    // 2. TG-NDA is needed for future research.
  }
  return equ_ptrs;
}

template <int dim>
std::unique_ptr<EigenBase<dim>> BuildEigenItr (
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr) {
  const std::string eigen_name(prm.get("eigen solver name"));
  std::unordered_map<std::string, int> mp = {{"pi", 0}};
  std::unique_ptr<EigenBase<dim>> eig;
  switch (mp[eigen_name]) {
    case 0:
      eig = std::unique_ptr<EigenBase<dim>>(
          new PowerIteration<dim>(prm, dat_ptr));
      break;
    default:
      break;
  }
  return eig;
}

template <int dim>
std::unique_ptr<MGBase<dim>> BuildMGItr (
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr) {
  const std::string mg_name(prm.get("mg solver name"));
  std::unordered_map<std::string, int> mp = {{"gs", 0}};
  std::unique_ptr<MGBase<dim>> mg;
  switch (mp[mg_name]) {
    case 0:
      mg = std::unique_ptr<MGBase<dim>>(new GaussSeidel<dim>(prm, dat_ptr));
      break;
    default:
      break;
  }
  return mg;
}

template <int dim>
std::unique_ptr<IGBase<dim>> BuildIGItr (
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr) {
  const std::string ig_name(prm.get("in group solver name"));
  std::unordered_map<std::string, int> mp = {{"si", 0}};
  std::unique_ptr<IGBase<dim>> ig;
  switch (mp[ig_name]) {
    case 0:
      ig = std::unique_ptr<IGBase<dim>>(new SourceIteration<dim>(prm, dat_ptr));
      break;
    default:
      break;
  }
  return ig;
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


// explicit instantiation for BuildAQ
template void bbuilders::BuildAQ<1> (const dealii::ParameterHandler&,
    std::unique_ptr<AQBase<1>> &);
template void bbuilders::BuildAQ<2> (const dealii::ParameterHandler&,
    std::unique_ptr<AQBase<2>> &);
template void bbuilders::BuildAQ<3> (const dealii::ParameterHandler&,
    std::unique_ptr<AQBase<3>> &);

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

template std::unique_ptr<IGBase<1>> bbuilders::BuildIGItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<1>> &);
template std::unique_ptr<IGBase<2>> bbuilders::BuildIGItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<2>> &);
template std::unique_ptr<IGBase<3>> bbuilders::BuildIGItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<3>> &);

template std::unique_ptr<MGBase<1>> bbuilders::BuildMGItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<1>> &);
template std::unique_ptr<MGBase<2>> bbuilders::BuildMGItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<2>> &);
template std::unique_ptr<MGBase<3>> bbuilders::BuildMGItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<3>> &);

template std::unique_ptr<EigenBase<1>> bbuilders::BuildEigenItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<1>> &);
template std::unique_ptr<EigenBase<2>> bbuilders::BuildEigenItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<2>> &);
template std::unique_ptr<EigenBase<3>> bbuilders::BuildEigenItr (
    const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<3>> &);

template std::unordered_map<std::string, std::unique_ptr<EquationBase<1>>>
    bbuilders::BuildEqu(const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<1>> &);
template std::unordered_map<std::string, std::unique_ptr<EquationBase<2>>>
    bbuilders::BuildEqu(const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<2>> &);
template std::unordered_map<std::string, std::unique_ptr<EquationBase<3>>>
    bbuilders::BuildEqu(const dealii::ParameterHandler &,
    std::shared_ptr<FundamentalData<3>> &);
