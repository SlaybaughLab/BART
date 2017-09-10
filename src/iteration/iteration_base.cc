#include "iteration_base.h"

template <int dim>
IterationBase<dim>::IterationBase (ParameterHandler &prm)
:
n_group(prm.get_integer("number of groups")),
do_nda(prm.get_bool("do NDA"))
{
}

template <int dim>
IterationBase<dim>::~IterationBase ()
{
}

template <int dim>
void IterationBase<dim>::do_iterations
(std::vector<Vector<double> > &sflx_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs)
{
}

// *****************************************************************************
// The following section are four member functions used to estimate difference
// between two vectors. The difference estimation is based on L1 norm.
// *****************************************************************************

template <int dim>
double IterationBase<dim>::estimate_phi_diff
(std::vector<PETScWrappers::MPI::Vector*> &phis_newer,
 std::vector<PETScWrappers::MPI::Vector*> &phis_older)
{
  AssertThrow (phis_newer.size ()== phis_older.size (),
               ExcMessage ("n_groups for different phis should be identical"));
  double err = 0.0;
  for (unsigned int i=0; i<phis_newer.size (); ++i)
  {
    PETScWrappers::MPI::Vector dif = *phis_newer[i];
    dif -= *phis_older[i];
    err = std::max (err, dif.l1_norm () / phis_newer[i]->l1_norm ());
  }
  return err;
}

template <int dim>
double IterationBase<dim>::estimate_phi_diff
(PETScWrappers::MPI::Vector* phi_newer,
 PETScWrappers::MPI::Vector* phi_older)
{
  PETScWrappers::MPI::Vector dif = *phi_newer;
  dif -= *phi_older;
  return std::max (err, dif.l1_norm () / phi_newer->l1_norm ());
}

template <int dim>
double IterationBase<dim>::estimate_phi_diff
(std::vector<Vector<double> > &phis1,
 std::vector<Vector<double> > &phis2)
{
  AssertThrow (phis1.size()==phis2.size(),
               ExcMessage("lengths of both vectors of solutions should be the same"));
  double err = 0.0;
  for (unsigned int g=0; g<phis1.size(); ++g)
  {
    Vector<double> dif = phis1[g];
    dif -= phis2[g];
    double local_numerator = dif.l1_norm ();
    double numerator = Utilities::MPI::sum (local_numerator, MPI_COMM_WORLD);
    double local_denorminator = phis1[g].l1_norm ();
    double denorminator = Utilities::MPI::sum (local_denorminator, MPI_COMM_WORLD);
    err = std::max(err, numerator/denorminator);
  }
  return err;
}

template <int dim>
double IterationBase<dim>::estimate_phi_diff
(Vector<double> &phi1, Vector<double> &phi2)
{
  Vector<double> dif = phi1;
  dif -= phi2;
  double local_numerator = dif.l1_norm ();
  double numerator = Utilities::MPI::sum (local_numerator, MPI_COMM_WORLD);
  double local_denorminator = phi1.l1_norm ();
  double denorminator = Utilities::MPI::sum (local_denorminator, MPI_COMM_WORLD);
  return numerator / denorminator;
}

template class IterationBase<2>;
template class IterationBase<3>;
