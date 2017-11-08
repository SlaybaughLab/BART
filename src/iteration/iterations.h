#ifndef __iterations_h__
#define __iterations_h__

#include <iostream>
#include <utility>
#include <vector>

#include "eigen_base.h"
#include "mg_base.h"
#include "ig_base.h"
#include "../equation/equation_base.h"


using namespace dealii;

//! This class provide an interface between BartDriver and specific iteration schemes.
/*!
 This class implements an interface between BartDriver<dim> and iteration schemes.
 The main purpose is such that one would not need to develop different functions
 according to different problem types.
 */
template <int dim>
class Iterations
{
public:
  /*!
   Class constructor.
   
   \param prm const ParameterHandler object.
   */
  Iterations (const ParameterHandler &prm);
  
  //! Class destructor.
  ~Iterations ();
  
  /*!
   A generic function to solve the problem for both eigenvalue and fixed-source
   problems. Inside this function, solving methods will be specified according to
   Iterations<dim>::is_eigen_problem the boolean.
   
   \note In fixed-source problems, eig_ptr points to a null pointer and will do 
   nothing.
   
   \param equ_ptrs A vector of shared_ptr's of EquationBase objects.
   \param ig_ptr A shared_ptr of downcast IGBase object.
   \param mg_ptr A shared_ptr of downcast MGBase object.
   \param eig_ptr A shared_ptr of downcast EigenBase object.
   \return Void.
   */
  void solve_problems
  (std::vector<Vector<double> > &sflxes_proc,
   std::vector<std_cxx11::shared_ptr<EquationBase<dim> > > &equ_ptrs,
   std_cxx11::shared_ptr<IGBase<dim> > ig_ptr,
   std_cxx11::shared_ptr<MGBase<dim> > mg_ptr,
   std_cxx11::shared_ptr<EigenBase<dim> > eig_ptr);

private:
  const std::string transport_name;//!< HO equation name.
  
  double keff;//!< keff living in current class.
  bool is_eigen_problem;//!< Boolean to determine if the problem is eigenvalue problem.
};

#endif	// define  __iterations_h__
