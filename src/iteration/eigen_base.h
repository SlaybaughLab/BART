#ifndef __eigen_base_h__
#define __eigen_base_h__

#include "iteration_base.h"

template <int dim>
class EigenBase : public IterationBase
{
public:
  EigenBase (ParameterHandler &prm);
  virtual ~EigenBase ();
  
  virtual void do_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs);
  
  virtual void eigen_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs);
  
  virtual void update_prev_sflxes_fiss_src_keff
  (std::vector<Vector<double> >&sflxes_proc);
  
  void get_keff (double &k);

protected:
  double estimate_fiss_src (std::vector<Vector<double> > &sflxes_proc);
  double estimate_k ();
  double estimate_k_err ();

  const double err_k_tol;
  const double err_phi_tol;

  double keff;
  double keff_prev;
  double fiss_src;
  double fiss_src_prev;

  std::vector<Vector<double> > sflxes_proc_prev_eigen;
  
  std_cxx11::shared_ptr<MGBase<dim> > mg_ptr;
}

#endif //__eigen_base_h__
