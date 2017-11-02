#ifndef __iteration_base_h__
#define __iteration_base_h__

#include <deal.II/base/parameter_handler.h>

#include "../equation/equation_base.h"

using namespace dealii;

template <int dim>
class IterationBase
{
public:
  IterationBase (const ParameterHandler &prm);
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
  
  const unsigned int n_group;
  const bool is_eigen_problem;
  const bool do_nda;
  
  double total_calculation_time; /**< total time for calculations+assemblies*/
  unsigned int ct_ho_iters; /**< HO iteration counts*/
  unsigned int ct_nda_iters; /**< NDA iteration counts*/
  
  ConditionalOStream pcout;
};

#endif // __iteration_base_h__
