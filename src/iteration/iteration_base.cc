#include "iteration_base.h"

template <int dim>
IterationBase<dim>::IterationBase (const ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> dat_ptr)
    :
    mat_vec_(dat_ptr->mat_vec),
    n_group_(prm.get_integer("number of groups")),
    is_eigen_problem_(prm.get_bool("do eigenvalue calculations")),
    do_nda_(prm.get_bool("do NDA")),
    ho_equ_name_(prm.get("transport model")) {}

template <int dim>
IterationBase<dim>::~IterationBase () {}

// *****************************************************************************
// The following section are four member functions used to estimate difference
// between two vectors. The difference estimation is based on L1 norm.
// *****************************************************************************

template <int dim>
double IterationBase<dim>::EstimatePhiDiff (
    std::vector<PETScWrappers::MPI::Vector*> &phis1,
    std::vector<PETScWrappers::MPI::Vector*> &phis2) {
  AssertThrow (phis1.size ()== phis2.size (),
               ExcMessage ("n_groups for different phis should be identical"));
  double err = 0.0;
  for (int i=0; i<phis1.size (); ++i)
    err = std::max (err, EstimatePhiDiff(phis1[i],phis2[i]));
  return err;
}

template <int dim>
double IterationBase<dim>::EstimatePhiDiff (
    dealii::PETScWrappers::MPI::Vector* phi1,
    dealii::PETScWrappers::MPI::Vector* phi2) {
  dealii::PETScWrappers::MPI::Vector dif = *phi1;
  dif -= *phi2;
  return dif.l1_norm () / phi1->l1_norm ();
}

template <int dim>
double IterationBase<dim>::EstimatePhiDiff (
    dealii::Vector<double> &phi1, dealii::Vector<double> &phi2) {
  dealii::Vector<double> dif = phi1;
  dif -= phi2;
  double local_numerator = dif.l1_norm ();
  double numerator = Utilities::MPI::sum (local_numerator, MPI_COMM_WORLD);
  double local_denorminator = phi1.l1_norm ();
  double denorminator = Utilities::MPI::sum (local_denorminator, MPI_COMM_WORLD);
  return numerator / denorminator;
}

template <int dim>
double IterationBase<dim>::EstimatePhiDiff (
    std::map<std::pair<int,int>, dealii::Vector<double>> &phis1,
    std::map<std::pair<int,int>, dealii::Vector<double>> &phis2) {
  AssertThrow (phis1.size()==phis2.size(),
      ExcMessage("lengths of both vectors of solutions should be the same"));
  double err = 0.0;
  for (int g=0; g<phis1.size(); ++g)
    err = std::max(err, EstimatePhiDiff(phis1[{g,0}], phis2[{g,0}]));
  return err;
}

template class IterationBase<1>
template class IterationBase<2>;
template class IterationBase<3>;
