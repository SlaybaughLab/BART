#ifndef __even_parity__
#define __even_parity__

#include "transport_base.h"



template<int dim>
class EvenParity : public TransportBase<dim>
{
public:
  EvenParity (ParameterHandler &prm);
  ~EvenParity ();
  
  void pre_assemble_cell_matrices
  (std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  void integrate_cell_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   FullMatrix<double> &cell_matrix,
   unsigned int &i_dir,
   unsigned int &g,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  void integrate_boundary_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   FullMatrix<double> &cell_matrix,
   unsigned int &i_dir,
   unsigned int &g);
  
  void integrate_interface_bilinear_form
  (typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,
   unsigned int &i_dir,
   unsigned int &g,
   FullMatrix<double> &vp_up,
   FullMatrix<double> &vp_un,
   FullMatrix<double> &vn_up,
   FullMatrix<double> &vn_un);
}

#endif // __even_parity__
