#include "bart_builder.h"

#include <unordered_map>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>

template <int dim>
BARTBuilder<dim>::BARTBuilder (dealii::ParameterHandler &prm) {
  SetParams (prm);
}

template <int dim>
BARTBuilder<dim>::~BARTBuilder () {}

template <int dim>
void BARTBuilder<dim>::SetParams (dealii::ParameterHandler &prm) {
  // for FE builder
  do_nda_ = prm.get_bool ("do nda");
  p_order_ = prm.get_integer("finite element polynomial degree");
  ho_discretization_ = prm.get ("ho spatial discretization");
  nda_discretization_ = prm.get ("nda spatial discretization");


  // for AQ builder
  aq_name_ = prm.get ("angular quadrature name");
  aq_order_ = prm.get_integer ("angular quadrature order");
}

template <int dim>
void BARTBuilder<dim>::BuildFESpaces (
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

      default:
        AssertThrow (false,
            dealii::ExcMessage("Invalid NDA discretization name"));
        break;
    }
  }
}

template <int dim>
void BARTBuilder<dim>::BuildAQ (dealii::ParameterHandler &prm,
    std::unique_ptr<AQBase<dim>> aq_ptr) {
  if (dim==1) {
    // AQBase implements 1D quadrature
    aq_ptr = std::unique_ptr<AQBase<dim>> (new AQBase<dim>(prm));
  } else if (dim>1) {
    std::unordered_map<std::string, int> aq_ind = {{"lsgc", 0}};

    switch (aq_ind[aq_name_]) {
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

template class BARTBuilder<1>;
template class BARTBuilder<2>;
template class BARTBuilder<3>;
