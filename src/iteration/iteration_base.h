#ifndef __iteration_base_h__
#define __iteration_base_h__

#include <deal.II/base/parameter_handler.h>

#include "../equation/equation_base.h"

using namespace dealii;

//! This class provides common functionalities for iteration related classes.
/*!
 \author Weixiong
 \date 2017/08~10
 \todo Implement iteration counts for every type of iteration in calculations.
 */
template <int dim>
class IterationBase
{
public:
  /*!
   Class constructor.
   
   \param prm Const dealii::ParameterHandler object.
   */
  IterationBase (const ParameterHandler &prm);
  
  //! Virtual class destructor.
  virtual ~IterationBase ();

protected:
  /** \brief Function to measure the relative difference between two sets of PETSc
   * Vectors
   *
   * \param Two sets of PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double estimate_phi_diff
  (std::vector<PETScWrappers::MPI::Vector*> &phis_newer,
   std::vector<PETScWrappers::MPI::Vector*> &phis_older);
  
  /** Function to measure the relative difference between two PETSc
   * Vectors
   *
   * \param Two PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double estimate_phi_diff
  (PETScWrappers::MPI::Vector* phi_newer,
   PETScWrappers::MPI::Vector* phi_older);
  
  /**
   * Function to measure the relative difference between two PETSc MPI Vectors.
   *
   * \param Two sets of PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double estimate_phi_diff
  (std::vector<Vector<double> > &phis_newer,
   std::vector<Vector<double> > &phis_older);
  
  /**
   * Function to measure the relative difference between two sets of deal.II 
   * Vectors. Note that though all vectors live on current processor, broadcast
   * is performed to obtain global results.
   *
   * \param Two deal.II Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double estimate_phi_diff
  (Vector<double> &phi_newer, Vector<double> &phi_older);
  
  const unsigned int n_group;//!< Number of groups.
  const bool is_eigen_problem;//!< Boolean to determine if it's eigenvalue problem.
  const bool do_nda;//!< Boolean to determine if NDA is used.
  
  double total_calculation_time; /**< total time for calculations+assemblies*/
  unsigned int ct_ho_iters; /**< HO iteration counts*/
  unsigned int ct_nda_iters; /**< NDA iteration counts*/
  
  //! ostream on processor with rank to be 0.
  /*!
   Details can be found <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/cl
   assConditionalOStream.html" style="color:blue"><b>here</b></a>.
   */
  ConditionalOStream pcout;
};

#endif // __iteration_base_h__
