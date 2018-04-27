#include "bart_builder.h"

#include <unordered_map>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>

namespace bbuilders {
template <int dim>
void BuildFESpaces (const dealii::ParameterHandler &prm,
    std::vector<dealii::FiniteElement<dim, dim>*> &fe_ptrs) {
  // getting parameter values
  const bool do_nda = prm.get_bool ("do nda");
  const int p_order = prm.get_integer("finite element polynomial degree");
  const std::string ho_discretization = prm.get ("ho spatial discretization");
  const std::string nda_discretization = prm.get ("nda spatial discretization");

  fe_ptrs.resize (do_nda ? 2 : 1);

  std::unordered_map<std::string, unsigned int> discretization_ind = {
      {"cfem", 0}, {"dfem", 1}, {"cmfd", 2}, {"rtk", 3}};

  switch (discretization_ind[ho_discretization]) {
    case 0:
      fe_ptrs.front () = new dealii::FE_Q<dim> (p_order);
      break;

    case 1:
      fe_ptrs.front () = new dealii::FE_DGQ<dim> (p_order);
      break;

    default:
      AssertThrow (false,
          dealii::ExcMessage("Invalid HO discretization name"));
      break;
  }

  if (do_nda) {
    switch (discretization_ind[nda_discretization]) {
      case 0:
        fe_ptrs.back () = new dealii::FE_Q<dim> (p_order);
        break;

      case 1:
        fe_ptrs.back () = new dealii::FE_DGQ<dim> (p_order);
        break;

      case 2:
        fe_ptrs.back () = new dealii::FE_DGQ<dim> (0);
        break;

      case 3:
        fe_ptrs.back () = new dealii::FE_RaviartThomas<dim> (p_order);
        break;

      default:
        AssertThrow (false,
            dealii::ExcMessage("Invalid NDA discretization name"));
        break;
    }
  }
}

template <int dim>
void BuildAQ (const dealii::ParameterHandler &prm,
    std::unique_ptr<AQBase<dim>> &aq_ptr) {
  // getting parameter values
  const std::string aq_name = prm.get ("angular quadrature name");
  if (dim==1) {
    // AQBase implements 1D quadrature
    aq_ptr = std::unique_ptr<AQBase<dim>> (new AQBase<dim>(prm));
  } else if (dim>1) {
    std::unordered_map<std::string, int> aq_ind = {{"lsgc", 0}};

    switch (aq_ind[aq_name]) {
      case 0:
        aq_ptr = std::unique_ptr<AQBase<dim>> (new LSGC<dim>(prm));
        break;

      default:
        AssertThrow (false,
            dealii::ExcMessage("Proper name is not given for AQ"));
        break;
    }
  }
}
}

// explicitly instantiate all builders using templates
// explicit instantiation for BuildFESpaces
template void bbuilders::BuildFESpaces<1> (const dealii::ParameterHandler&,
    std::vector<dealii::FiniteElement<1, 1>*> &);
template void bbuilders::BuildFESpaces<2> (const dealii::ParameterHandler&,
    std::vector<dealii::FiniteElement<2, 2>*> &);
template void bbuilders::BuildFESpaces<3> (const dealii::ParameterHandler&,
    std::vector<dealii::FiniteElement<3, 3>*> &);

// explicit instantiation for BuildAQ
template void bbuilders::BuildAQ<1> (const dealii::ParameterHandler&,
    std::unique_ptr<AQBase<1>> &);
template void bbuilders::BuildAQ<2> (const dealii::ParameterHandler&,
    std::unique_ptr<AQBase<2>> &);
template void bbuilders::BuildAQ<3> (const dealii::ParameterHandler&,
    std::unique_ptr<AQBase<3>> &);
