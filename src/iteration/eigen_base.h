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
  
  /*!
   Function to update \f$\phi\f$, fission source and \f$k_\mathrm{eff}\f$ from
   previous iteration with current values.
   
   \param sflxes_proc \f$\phi\f$ for all groups on current processor.
   \return Void
   */
  virtual void update_prev_sflxes_fiss_src_keff
  (std::vector<Vector<double> >&sflxes_proc);
  
  /*!
   Function to get \f$k_\mathrm{eff}\f$ value from eigenvalue iterations.
   
   \return \f$k_\mathrm{eff}\f$ value (double type).
   */
  double get_keff ();

protected:
  /*!
   Function to initialize \f$k_\mathrm{eff}\f$ and fission source. 
   \f$k_\mathrm{eff}\f$ is initilized with unit value while fission source is
   initialize with unit-value scalar fluxes for all groups.
   
   \param sflxes_proc A vector of scalar fluxes of all groups living on current
   processor.
   \param equ_ptr A pointer of EquationBase object.
   \return Void.
   */
  void initialize_fiss_process
  (std::vector<Vector<double> > &sflxes_proc,
   std_cxx11::shared_ptr<EquationBase<dim> > &equ_ptr);
  
  /*!
   Function to estimate fission source and \f$k_\mathrm{eff}\f$ values.
   
   \param sflxes_proc
   \return Void.
   */
  void calculate_fiss_src_keff
  (std::vector<Vector<double> > &sflxes_proc,
   std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr);
  
  /*!
   Function to calculate \f$k_\mathrm{eff}\f$.
   
   \return \f$k_\mathrm{eff}\f$ value.
   */
  double estimate_k ();
  
  /*!
   Function to calculate difference between \f$k_\mathrm{eff}\f$ and 
   \f$k_\mathrm{eff,\ prev}\f$.
   
   \return The relative difference (double type).
   */
  double estimate_k_diff ();

  const double err_k_tol;//!< Tolerance for convergence check on \f$k_\mathrm{eff}\f$
  const double err_phi_tol;//!< Tolerance for convergence check on \f$\phi\f$

  double keff;//!< \f$k_\mathrm{eff}\f$ value.
  double keff_prev;//!< \f$k_\mathrm{eff}\f$ value from previous eigen iteration.
  double fiss_src;//!< Fission source.
  double fiss_src_prev;//!< Fission source from previous eigen iteration.

  //! \f$\phi\f$ from previous eigen iteration for all groups
  std::vector<Vector<double> > sflxes_proc_prev_eigen;
};

#endif //__eigen_base_h__
