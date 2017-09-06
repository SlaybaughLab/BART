#ifndef __in_group_base_h__
#define __in_group_base_h__

template <int dim>
class InGroupBase : public IterationBase<dim>
{
public:
  InGroupBase ();
  virtual ~InGroupBase();
  
  // has to be provided
  virtual void solve_in_group
  (std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells,
   std::vector<bool> &is_cell_at
   std_cxx11::shared_ptr<EquationBase<dim> > equ_ptrs,
   unsigned int &g);
  
protected:
  const double err_phi_tol;
  
  std::vector<Vector<double> > sflx_proc_old;
};

template <int dim>
class SourceIteration : public InGroupBase<dim>
{
public:
  SourceIteration ();
  ~SourceIteration ();
  
  void solve_in_group
  (std_cxx11::shared_ptr<EquationBase<dim> > equ_ptrs,
   unsigned int &g);
};

#endif //__in_group_base_h__
