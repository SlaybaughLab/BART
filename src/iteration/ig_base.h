#ifndef __ig_base_h__
#define __ig_base_h__

#include "iteration_base.h"

using namespace dealii;

//! This class provides in-group solving scheme.
/*!
 This class serves as the base class for in group solvers. The main (only)
 functionality is to provide an abstract function for solve in a specific group.
 Overriding has to be provided in derived class.
 
 \author Weixiong Zheng
 \date 2017/09
 */
template <int dim>
class IGBase : public IterationBase<dim>
{
public:
  /*!
   Class constructor.
   
   \param prm Const ParameterHandler object.
   */
  IGBase (const ParameterHandler &prm);
  
  //! Class destructor.
  virtual ~IGBase();
  
  /*!
   Abstract function for in group solving functionality. For instance, one could
   provid overriding as one-pass solve for diffusion/NDA, or iteration based in
   group solving such as source iteration.
   
   \param sflxes_proc A vector of all scalar fluxes on current processor.
   \param equ_ptr A pointer to the equation of interest.
   \param g Group index.
   */
  virtual void solve_in_group
  (std::vector<Vector<double> > &sflxes_proc,
   std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
   unsigned int &g);
  
protected:
  const double err_phi_tol;//!< Tolerance used in iterative in group scheme for convergence check.
  
  /*!
   Scalar flux value from provious generation.
   
   \todo Needs modification if anisotropic scattering is to be implemented.
   Specifically, a quantities like std::vector<Vector<double> > has to be provided
   to contain all moments.
   */
  Vector<double> sflx_proc_prev_ig;
};

//! This class provides source iteration scheme for in group solve.
/*!
 \author Weixiong Zheng
 \date 2017/08
 */
template <int dim>
class SourceIteration : public IGBase<dim>
{
public:
  /*!
   Class constructor.
   
   \param prm Const ParameterHandler object.
   */
  SourceIteration (const ParameterHandler &prm);
  
  //! Class destructor.
  ~SourceIteration ();
  
  /*!
   A function to solve SN in group using source iteration scheme.
   
   \param sflxes_proc Scalar fluxes for all groups living on current processor.
   \param equ_ptr A pointer to the equation of interest.
   \param g Group index.
   */
  void solve_in_group
  (std::vector<Vector<double> > &sflxes_proc,
   std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
   unsigned int &g);
};

#endif //__ig_base_h__
