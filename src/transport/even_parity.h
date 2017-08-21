#ifndef __even_parity__
#define __even_parity__

#include "equation_base.h"

template<int dim>
class EvenParity : public EquationBase<dim>
{
public:
  EvenParity (ParameterHandler &prm,
              const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
              const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
              const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  ~EvenParity ();
  
  void pre_assemble_cell_matrices
  (const std_cxx11::shared_ptr<FEValues<dim> > fv,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  void integrate_cell_bilinear_form
  (const std_cxx11::shared_ptr<FEValues<dim> > fv,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   FullMatrix<double> &cell_matrix,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  void integrate_boundary_bilinear_form
  (const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,/*face number*/
   FullMatrix<double> &cell_matrix,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  void integrate_interface_bilinear_form
  (const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
   const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf_nei,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   typename DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
   unsigned int &fn,/*concerning face number in local cell*/
   FullMatrix<double> &vp_up,
   FullMatrix<double> &vp_un,
   FullMatrix<double> &vn_up,
   FullMatrix<double> &vn_un,
   const unsigned int &g,
   const unsigned int &i_dir);
  
  void generate_ho_fixed_source
  (std::vector<PETScWrappers::MPI::Vector*> &vec_ho_fixed_rhs,
   std::vector<Vector<double> > &sflx_this_proc);
  
  void generate_ho_rhs
  (std::vector<PETScWrappers::MPI::Vector*> &vec_ho_rhs,
   std::vector<PETScWrappers::MPI::Vector*> &vec_ho_fixed_rhs,
   std::vector<Vector> &sflx_this_proc);
  
private:
  double c_penalty;
  std::vector<double> tensor_norms;
};

#endif // __even_parity__
