#include "mg_base.h"

template <int dim>
MGBase<dim>::MGBase (const ParameterHandler &prm)
:
IterationBase<dim> (prm),
g_thermal(prm.get_integer("thermal group boundary")),
err_phi_tol(1.0e-5)
{
  AssertThrow (g_thermal<this->n_group && g_thermal>=0,
               ExcMessage("Invalid thermal upper boundary"));
  
  // TODO: needs change if anisotropic scattering is needed
  sflxes_proc_prev_mg.resize (this->n_group);
}

template <int dim>
MGBase<dim>::~MGBase ()
{
}

template <int dim>
void MGBase<dim>::do_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
 std_cxx11::shared_ptr<IGBase<dim> > ig_ptr)
{
  // this function will be called only when it's not an eigenvalue problem. Copy
  // the following assertion in ALL derived class of MGBase.
  AssertThrow (!this->is_eigen_problem,
               ExcMessage("This function shall not be called if it's eigen prob."));
  
  // override this function per derived class
}

template <int dim>
void MGBase<dim>::mg_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
 std_cxx11::shared_ptr<IGBase<dim> > ig_ptr)
{// this function needs to be overridden if JFNK is desired
  if (this->do_nda)
    AssertThrow (equ_ptrs.back()->get_equ_name()=="nda" &&
                 equ_ptrs.front()->get_equ_name()!="nda",
                 ExcMessage("Check equation names"));
  // solve for epithermal and fast groups
  nonthermal_solves (sflxes_proc, equ_ptrs, ig_ptr);
  // thermal group iterations
  thermal_iterations (sflxes_proc, equ_ptrs, ig_ptr);
}

template <int dim>
void MGBase<dim>::nonthermal_solves
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
 std_cxx11::shared_ptr<IGBase<dim> > ig_ptr)
{
}

template <int dim>
void MGBase<dim>::thermal_iterations
(std::vector<Vector<double> > &sflxes_proc,
 std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
 std_cxx11::shared_ptr<IGBase<dim> > ig_ptr)
{
}

template class MGBase<2>;
template class MGBase<3>;
