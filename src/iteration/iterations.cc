#include "iterations.h"
#include "../common/bart_builder.h"

template <int dim>
Iterations<dim>::Iterations (const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    is_eigen_problem_(prm.get_bool("do eigenvalue calculations")) {
  if (is_eigen_problem_) {
    eig_ptr_ = bbuilders::BuildEigenItr (prm, dat_ptr);
  } else {
    mg_ptr_ = bbuilders::BuildMGItr (prm, dat_ptr);
  }
}

template <int dim>
Iterations<dim>::~Iterations () {}

template <int dim>
void Iterations<dim>::DoIterations (std::unordered_map<std::string,
    std::unique_ptr<EquationBase<dim>>> &equ_ptrs) {
  if (is_eigen_problem_) {
    eig_ptr_->DoIterations (equ_ptrs);
    keff_ = eig_ptr_->GetKeff ();
  } else {
    mg_ptr_->DoIterations (equ_ptrs);
  }
}

template <int dim>
double Iterations<dim>::GetKeff () const {
  AssertThrow (is_eigen_problem_,
      dealii::ExcMessage("Problem is not eigenvalue problem"));
  return keff_;
}

// explicit instantiation to avoid linking error
template class Iterations<1>;
template class Iterations<2>;
template class Iterations<3>;
