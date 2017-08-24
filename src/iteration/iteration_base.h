#ifndef __iteration_base_h__
#define __iteration_base_h__

using namespace dealii;

template <int dim>
class IterationBase
{
public:
  IterationBase ();
  virtual ~IterationBase ();
  
  virtual void do_iterations ();
  virtual void assemble_system
  (std::vector<PETScWrappers::MPI::SparseMatrix*> sys_mats);
  
protected:
  /** \brief Function to measure the relative difference between two sets of PETSc
   * Vectors
   *
   * \parameters Two sets of PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double estimate_phi_diff
  (std::vector<PETScWrappers::MPI::Vector*> &phis_newer,
   std::vector<PETScWrappers::MPI::Vector*> &phis_older);
  
  /** \brief Function to measure the relative difference between two PETSc
   * Vectors
   *
   * \parameters Two PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double estimate_phi_diff
  (PETScWrappers::MPI::Vector* phi_newer,
   PETScWrappers::MPI::Vector* phi_older);
  
  /** \brief Function to measure the relative difference between two sets of PETSc
   * Vectors
   *
   * \parameters Two sets of PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double estimate_phi_diff
  (std::vector<Vector<double> > &phis_newer,
   std::vector<Vector<double> > &phis_older);
  
  /** \brief Function to measure the relative difference between two PETSc
   * Vectors
   *
   * \parameters Two PETSc Vectors with the same length
   * \return relative difference in vector l1 norm measure
   */
  double estimate_phi_diff
  (Vector<double> &phi_newer, Vector<double> &phi_older);
  
  std_cxx11::shared_ptr<EquationBase<dim> > trm_ptr;
  std_cxx11::shared_ptr<EquationBase<dim> > nda_ptr;
  
  void initialize_equations ();
  
  std_cxx11::shared_ptr<EquationBase<dim> > tra_ptr;
  std_cxx11::shared_ptr<EquationBase<dim> > nda_ptr;
  
  double total_calculation_time; /**< total time for calculations including assembly of rhs*/
  unsigned int ct_ho_iters; /**< HO iteration counts*/
  unsigned int ct_nda_iters; /**< NDA iteration counts*/
  
  std::vector<Vector<double> > sflx_proc;
  std::vector<Vector<double> > sflx_proc_prev_gen;
  std::vector<Vector<double> > lo_sflx_proc;
}

#endif // __iteration_base_h__
