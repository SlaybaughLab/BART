#ifndef __eigen_base_h__
#define __eigen_base_h__

#include "iteration_base.h"
#include "mg_base.h"

using namespace dealii;

//! This class provides eigenvalue calculations foundations.
/*!
 This class is the base class of eigenvalue calculations. It will operate the
 iteration at a high level without touching specifics in the equations. Alternatively,
 it interfaces with MGBase objects iteratively.
 
 \author Weixiong Zheng
 \date 2017/09
 */
template <int dim>
class EigenBase : public IterationBase<dim>
{
public:
  EigenBase (const ParameterHandler &prm);
  virtual ~EigenBase ();
  
  /*!
   This function calls eigen_iterations to performs eigenvalue iterations.
   
   \todo Add NDA functionality.
   */
  virtual void do_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
   std_cxx11::shared_ptr<MGBase<dim> > mg_ptr);
  
  /*!
   This virtual function performs eigenvalue iteration after overriding with the
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > >::back().
   
   \return Void.
   */
  virtual void eigen_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
   std_cxx11::shared_ptr<MGBase<dim> > mg_ptr);
  
  virtual void update_prev_sflxes_fiss_src_keff
  (std::vector<Vector<double> >&sflxes_proc);
  
  void get_keff (double &k);

protected:
  void initialize_fiss_process
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs);
  
  void calculate_fiss_src_keff
  (std::vector<Vector<double> > &sflxes_proc,
   std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr);
  
  double estimate_fiss_src (std::vector<Vector<double> > &sflxes_proc);
  double estimate_k ();
  double estimate_k_diff ();

  const double err_k_tol;
  const double err_phi_tol;

  double keff;
  double keff_prev;
  double fiss_src;
  double fiss_src_prev;

  std::vector<Vector<double> > sflxes_proc_prev_eigen;
};

#endif //__eigen_base_h__
