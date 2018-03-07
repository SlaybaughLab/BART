#include "bart_builder.h"

template <int dim>
BartBuilder<dim>::BartBuilder (dealii::ParameterHandler &prm) {
  InitCoefs (prm);
}

template <int dim>
BartBuilder<dim>::~BartBuilder () {
}

BartBuilder<dim>::InitCoefs (dealii::ParameterHandler &prm) {
  do_nda = prm.get_bool ("do NDA");
  p_order = prm.get_integer("finite element polynomial degree");
  ho_discretization = prm.get ("ho spatial discretization");
  nda_discretization = prm.get ("nda spatial discretization");
}

template <int dim>
BartBuilder<dim>::BuildFESpaces (
    std::std::vector<FiniteElement<dim, dim>*> &fe_ptrs) {
  fe_ptrs.resize (do_nda ? 2 : 1);
  switch (ho_discretization) {
    case "dfem" :
      fe_ptrs.front () = new FE_DGQ<dim> (p_order);
      break;

    // by default, it's "cfem"
    default:
      fe_ptrs.front () = new FE_Q<dim> (p_order);
      break;
  }

  if (do_nda) {
    switch (nda_discretization) {
      case "dfem":
        fe_ptrs.back () = new FE_DGQ<dim> (p_order);
        break;

      case "cmfd":
        fe_ptrs.back () = new FE_DGQ<dim> (0);
        break;

        // by default, it's "cfem"
      default:
        fe_ptrs.back () = new FE_Q<dim> (p_order);
        break;
    }
  }
}

template class BartBuilder<1>;
template class BartBuilder<2>;
template class BartBuilder<3>;
