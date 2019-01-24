#include "diffusion.h"

template <int dim>
Diffusion<dim>::Diffusion(const std::string &equation_name,
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    EquationBase<dim>(equation_name, prm, dat_ptr) {}

template <int dim>
Diffusion<dim>::~Diffusion () {}

template <int dim>
void Diffusion<dim>::PreassembleCellMatrices() {
  fv_->reinit(dat_ptr_->local_cells[0]);

  for (int qi = 0; qi < n_q_; ++qi) {
    pre_collision_[qi] = dealii::FullMatrix<double>(dofs_per_cell_,
                                                    dofs_per_cell_);
    pre_streaming_[qi] = dealii::FullMatrix<double>(dofs_per_cell_,
                                                    dofs_per_cell_);
    for (int i = 0; i < dofs_per_cell_; ++i) {
      for (int j = 0; j < dofs_per_cell_; ++j) {
        pre_collision_[qi](i,j) =
            fv_->shape_value(i,qi) * fv_->shape_value(j,qi);
        pre_streaming_[qi](i,j) =
            fv_->shape_grad(i,qi) * fv_->shape_grad(j,qi);
      }
    }
  }
}

template <int dim>
void Diffusion<dim>::GenerateMoments (
    std::map<std::tuple<int,int,int>, dealii::Vector<double>> &moments,
    std::map<std::tuple<int,int,int>, dealii::Vector<double>> &moments_prev,
    const int &g) {
  auto key = std::make_tuple(g,0,0);
  moments_prev[key] = moments[key];
  // generate moments
  moments[key] = *mat_vec_->sys_flxes[equ_name_][0];
}

template <int dim>
void Diffusion<dim>::IntegrateCellBilinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::FullMatrix<double> &cell_matrix,
    const int &g,
    const int &) {
  int material_id = cell->material_id ();
  auto sigma_r = xsec_->sigt.at(material_id)[g] - xsec_->sigs.at(material_id)(g,g);
  auto diffusion_coef = xsec_->diffusion_coef.at(material_id)[g];
  
  for (int q = 0; q < n_q_; ++q) {
    for (int i = 0; i < dofs_per_cell_; ++i) {
      for (int j = 0; j < dofs_per_cell_; ++j) {
        cell_matrix(i,j) += (pre_streaming_[q](i,j) * diffusion_coef                             
                             +
                             pre_collision_[q](i,j) * sigma_r) * fv_->JxW(q);
      }
    }
  }
}

template <int dim>
void Diffusion<dim>::IntegrateScatteringLinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &) {
  int material_id = cell->material_id ();

  // Scalar flux times sigma_s, at each quadrature point
  std::vector<double> cell_scatter_flux(n_q_);

  // Iterate over groups to populate cell_scatter_flux
  for (int group_in = 0; group_in < n_group_; ++group_in) {
    std::vector<double> group_cell_scatter_flux (n_q_);

    // Get needed cross-section
    auto sigma_s = xsec_->sigs.at(material_id)(group_in,g);
    
    // Get group scalar flux and populate group_cell_scatter_flux
    fv_->get_function_values(
        mat_vec_->moments[equ_name_][std::make_tuple(group_in, 0, 0)],
        group_cell_scatter_flux);

    // Sum over all groups at each quadrature point if sigma_s is non-zero
    if (g != group_in && sigma_s > bconst::kSmall) {
      for (int q = 0; q < n_q_; ++q)
        cell_scatter_flux[q] += sigma_s * group_cell_scatter_flux[q];
    }
  }

  // Integrate by summing over quadrature points (multiplied by Jacobian) and
  // add to the appropriate element of the RHS vector
  for (int q = 0; q < n_q_; ++q) {
    cell_scatter_flux[q] *= fv_->JxW(q);
    for (int i = 0; i < dofs_per_cell_; ++i)
      cell_rhs(i) += fv_->shape_value(i,q) * cell_scatter_flux[q];
  }
}

template <int dim>
void Diffusion<dim>::IntegrateCellFixedLinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &) {
  int material_id = cell->material_id ();

  // Fixed source value at each quadrature point in the cell
  std::vector<double> cell_q(n_q_);

  if (!is_eigen_problem_) {

    // Q value at each quadrature point is given by the material properties
    cell_q = std::vector<double>(n_q_, xsec_->q.at(material_id)[g]);
    
  } else if (xsec_->is_material_fissile.at(material_id)) {
    
    for (int group_in = 0; group_in < n_group_; ++group_in) {
      std::vector<double> group_cell_scalar_flux (n_q_);
      
      // Get group cell scalar flux
      fv_->get_function_values (
          mat_vec_->moments[equ_name_][std::make_tuple(group_in, 0, 0)],
          group_cell_scalar_flux);

      auto scaled_fission_transfer =
          scaled_fiss_transfer_.at(material_id)(group_in, g);
      
      for (int q = 0; q < n_q_; ++q)
        cell_q[q] += scaled_fission_transfer * group_cell_scalar_flux[q];
    }
  }
  
  // Integrate and add to RHS
  for (int q=0; q<n_q_; ++q) {
    
    cell_q[q] *= fv_->JxW(q);
    
    for (int i = 0; i < dofs_per_cell_; ++i)
      cell_rhs(i) += fv_->shape_value(i,q) * cell_q[q];
  }
}

template <int dim>
void Diffusion<dim>::IntegrateBoundaryBilinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &,
    const int &,/*face number*/
    dealii::FullMatrix<double> &cell_matrix,
    const int &,
    const int &) {

  // Vacuum boundary conditions
  if (!have_reflective_bc_) {
    for (int q = 0; q < n_qf_; ++q) {      
      for (int i = 0; i < dofs_per_cell_; ++i) {
        for (int j = 0; j < dofs_per_cell_; ++j) {
          cell_matrix(i,j) += (fvf_->shape_value(i,q) *
                               fvf_->shape_value(j,q) *
                               0.25 * fvf_->JxW(q));
        }
      }
    }
  }
}

template class Diffusion<1>;
template class Diffusion<2>;
template class Diffusion<3>;
