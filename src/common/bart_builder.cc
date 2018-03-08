#include "bart_builder.h"

#include <unordered_map>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>

template <int dim>
BartBuilder<dim>::BartBuilder (dealii::ParameterHandler &prm) {
  SetParams (prm);
}

template <int dim>
BartBuilder<dim>::~BartBuilder () {
}

template <int dim>
void BartBuilder<dim>::SetParams (dealii::ParameterHandler &prm) {
  do_nda_ = prm.get_bool ("do nda");
  p_order_ = prm.get_integer("finite element polynomial degree");
  ho_discretization_ = prm.get ("ho spatial discretization");
  nda_discretization_ = prm.get ("nda spatial discretization");
}

template <int dim>
void BartBuilder<dim>::BuildFESpaces (
    std::vector<dealii::FiniteElement<dim, dim>*> &fe_ptrs) {
  fe_ptrs.resize (do_nda_ ? 2 : 1);

  std::unordered_map<std::string, unsigned int> discretization_ind = {
      {"cfem", 0}, {"dfem", 1}, {"cmfd", 2}, {"rtk", 3}};

  switch (discretization_ind[ho_discretization_]) {
    case 0:
      fe_ptrs.front () = new dealii::FE_Q<dim> (p_order_);
      break;

    case 1:
      fe_ptrs.front () = new dealii::FE_DGQ<dim> (p_order_);
      break;

    // by default, it's "cfem"
    default:
      AssertThrow (false,
          dealii::ExcMessage("Invalid HO discretization name"));
      break;
  }

  if (do_nda_) {
    switch (discretization_ind[nda_discretization_]) {
      case 0:
        fe_ptrs.back () = new dealii::FE_Q<dim> (p_order_);
        break;

      case 1:
        fe_ptrs.back () = new dealii::FE_DGQ<dim> (p_order_);
        break;

      case 2:
        fe_ptrs.back () = new dealii::FE_DGQ<dim> (0);
        break;

      case 3:
        fe_ptrs.back () = new dealii::FE_RaviartThomas<dim> (p_order_);
        break;

      // by default, it's "cfem"
      default:
        AssertThrow (false,
            dealii::ExcMessage("Invalid NDA discretization name"));
        break;
    }
  }
}

template class BartBuilder<1>;
template class BartBuilder<2>;
template class BartBuilder<3>;
