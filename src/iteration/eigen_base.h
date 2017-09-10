#ifndef __eigen_base_h__
#define __eigen_base_h__

#include "iteration_base.h"

template <int dim>
class EigenBase : public IterationBase
{
public:
  EigenBase (ParameterHandler &prm);
  virtual ~EigenBase ();
  
  virtual void do_iterations ();
  
  virtual void eigen_iterations ();
  
  double get_keff ();
  
protected:
  double estimate_fission_source (std::vector<Vector<double> > &sflx_proc);
  double estimate_k (double &fiss_src, double &fiss_src_prev, double &k_prev);
  double estimate_k_err (double &k, double &k_prev);
  
  const double err_k_tol;
  const double err_phi_tol;
  
  std_cxx11::shared_ptr<MGBase<dim> > mg_ptr;
  
  double keff;
  double keff_prev;
  double fiss_source;
  double fiss_source_prev;
  
  std::vector<Vector<double> > sflxes_proc_prev_fiss;
}

#endif //__eigen_base_h__
