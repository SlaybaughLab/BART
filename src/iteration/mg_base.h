#ifndef __mg_base_h__
#define __mg_base_h__

using namespace dealii;

template <int dim>
class MGBase : public IterationBase<dim>
{
public:
  MGBase ();
  virtual ~MGBase ();
  
  virtual void mg_iterations ();
  virtual void generate_system_matrices
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats);
  virtual void generate_group_rhses
  (std::vector<PETScWrappers::MPI::Vector*> &group_rhses, unsigned int &g);
  virtual void iteration_over_groups
  (std::vector<PETScWrappers::MPI::Vector*> &group_rhses);
}

#endif//__mg_base_h__
