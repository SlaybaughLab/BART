#ifndef __gauss_seidel_h__
#define __gauss_seidel_h__

#include "mg_base.h"

//! This class provides Gauss-Seidel scheme for MG calculations.
/*!
 This class implements Gauss-Seidel scheme for MG calculations. Specifically,
 MGBase<dim>::nonthermal_solves and MGBase<dim>::thermal_iterations are overriden
 using Gauss-Seidel iterations.
 
 \author Weixiong Zheng
 \date 2017/09
 */
template <int dim>
class GaussSeidel : public MGBase<dim>
{
public:
  /*!
   Class constructor.
   
   \param prm const ParameterHandler object.
   */
  GaussSeidel (const ParameterHandler &prm);
  
  //! Class destructor.
  ~GaussSeidel ();
  
  /*!
   A function overriding MGBase<dim>::nonthermal_solves. The function implements
   Gauss-Seidel scheme with one-pass serial solves (energy groupwise) from highest
   energy group until reaching out to the last group before the thermal group 
   boundary.
   
   \param sflxes_proc Scalar fluxes living on current processor for all groups.
   \param equ_ptrs A vector of pointers of EquationBase<dim> object.
   \param ig_ptr A pointer of IGBase<dim> object.
   \return Void.
   */
  void nonthermal_solves
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
  
  /*!
   A function overriding MGBase<dim>::thermal_iterations. The function implements
   Gauss-Seidel scheme for solving transport problems in thermal groups 
   iteratively. The starting group is the highest energy group affected by 
   upscattering.
   
   \param sflxes_proc Scalar fluxes living on current processor for all groups.
   \param equ_ptrs A vector of pointers of EquationBase<dim> object.
   \param ig_ptr A pointer of IGBase<dim> object.
   \return Void.
   */
  void thermal_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
};

#endif //__gauss_Seidel_h__
