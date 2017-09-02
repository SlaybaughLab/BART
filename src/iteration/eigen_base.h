#ifndef __eigen_base_h__
#define __eigen_base_h__

template <int dim>
class EigenBase : public IterationBase
{
public:
  EigenBase ();
  virtual ~EigenBase ();
  
  virtual void do_iterations ();
  virtual void eigen_iterations ();
  
  double get_keff ();
  
protected:
  double estimate_fission_source (std::vector<Vector<double> > &sflx_proc);
  double estimate_k (double &fiss_src, double &fiss_src_prev_gen, double &k_prev);
  double estimate_k_err (double &k, double &k_prev);
  
  const double err_k_tol;
  const double err_phi_tol;
  
  double keff;
  double keff_prev;
  
  std::vector<Vector<double> > sflx_proc_old;
}

#endif //__eigen_base_h__
