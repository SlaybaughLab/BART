#ifndef __gauss_sidel_h__
#define __gauss_sidel_h__

#include "mg_base.h"

template <int dim>
class GaussSidel : public MGBase<dim>
{
public:
  GaussSidel (const ParameterHandler &prm);
  ~GaussSidel ();
  
  void nonthermal_solves
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
  
  void thermal_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
  
  void do_iterations
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr);
};

#endif //__gauss_sidel_h__
