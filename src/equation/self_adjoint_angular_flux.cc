#include "self_adjoint_angular_flux.h"
#include <deal.II/dofs/dof_accessor.h>

template <int dim>
SelfAdjointAngularFlux<dim>::SelfAdjointAngularFlux(
    const std::string equation_name,
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &data_ptr)
    : EquationBase<dim>(equation_name, prm, data_ptr) {}

/*
 * =============================================================================
 * PUBLIC FUNCTIONS
 * =============================================================================
 */

template <int dim>
void SelfAdjointAngularFlux<dim>::AssembleLinearForms (const int &g) {
  if (have_reflective_bc_) {
    for (int dir = 0; dir < n_dir_; ++dir) {
      global_angular_flux_[this->GetCompInd(g, dir)] =
          *(mat_vec_->sys_flxes[equ_name_][this->GetCompInd(g, dir)]);         
    }
  }
  EquationBase<dim>::AssembleLinearForms(g);
}

template<int dim>
void SelfAdjointAngularFlux<dim>::IntegrateBoundaryBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &,
      const int &,
      dealii::FullMatrix<double> &cell_matrix,
      const int &,
      const int &dir) {
  
  // Get the boundary ID and the normal vector
  const dealii::Tensor<1, dim> normal_vector = fvf_->normal_vector(0);

  double normal_dot_omega = normal_vector * omega_[dir];
  
  if (normal_dot_omega > 0) {
    // Integrate bilinear term
    for (int q = 0; q < n_qf_; ++q) {
      for (int i = 0; i < dofs_per_cell_; ++i) {
        for (int j = 0; j < dofs_per_cell_; ++j) {
          cell_matrix(i,j) += normal_dot_omega *
                              fvf_->shape_value(i, q) *
                              fvf_->shape_value(j, q) *
                              fvf_->JxW(q);
        }
      }
    }
  }
}

template<int dim>
void SelfAdjointAngularFlux<dim>::IntegrateBoundaryLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const int &fn,/*face number*/
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) {

  // Get boundary id to determine if we have reflective BCs
  int boundary_id = cell->face(fn)->boundary_id();

  if (have_reflective_bc_ && is_reflective_bc_.at(boundary_id)) {
  
    // Get normal vector and dot product with omega
    const dealii::Tensor<1, dim> normal_vector = fvf_->normal_vector(0);
    double normal_dot_omega = normal_vector * omega_[dir];

    if (normal_dot_omega < 0) {
      // Retrieve previous iteration angular flux and get local values for the
      // reflected angle
      int reflected_dir = this->ref_dir_ind_[std::make_pair(boundary_id, dir)];
      
      std::vector<double> angular_flux(n_qf_);
      fvf_->get_function_values(
          global_angular_flux_[this->GetCompInd(g, reflected_dir)],
          angular_flux);     
      
      for (int q = 0; q < n_qf_; ++q) {
        for (int i = 0; i < dofs_per_cell_; ++i) {
          cell_rhs(i) -= normal_dot_omega *
                         fvf_->shape_value(i, q) *
                         angular_flux[q] *
                         fvf_->JxW(q);
        }
      }
    }
  }
}

template<int dim>
void SelfAdjointAngularFlux<dim>::IntegrateCellBilinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::FullMatrix<double> &cell_matrix,
      const int &g,
      const int &dir) {
  // Get material id for the given cell and cross-sections
  int material_id = cell->material_id();
  auto sigma_t = xsec_->sigt.at(material_id)[g];
  auto inv_sigma_t = xsec_->inv_sigt.at(material_id)[g];

  // Integrate and add both bilinear terms using precomputed values
  for (int q = 0; q < n_q_; ++q) {
    for (int i = 0; i < dofs_per_cell_; ++i) {
      for (int j = 0; j < dofs_per_cell_; ++j) {

        cell_matrix(i, j) += (pre_streaming_[{dir, q}](i, j) * inv_sigma_t
                              +
                              pre_collision_[q](i, j) * sigma_t) * fv_->JxW(q);
      }
    }
  }
}

template<int dim>
void SelfAdjointAngularFlux<dim>::IntegrateCellFixedLinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &dir) {
  // get material id for the given cell
  int material_id = cell->material_id();
  
  // SAAF has two fixed source terms, the first is proportional to q/4pi, the
  // second contains an additional division by sigma_t. These two vectors hold
  // these values, respectively, for each quadrature point.
  std::vector<double> cell_q(n_q_);
  std::vector<double> cell_q_over_total(n_q_);
  
  if (!is_eigen_problem_) {
    
    // Fill the two vectors with the appropriate q/4pi and q/4pi*sigma_t values
    auto q_per_ster = xsec_->q_per_ster.at(material_id)[g];
    std::fill(cell_q.begin(), cell_q.end(),
              q_per_ster);
    std::fill(cell_q_over_total.begin(), cell_q_over_total.end(),
              q_per_ster * xsec_->inv_sigt.at(material_id)[g]);
    
  } else if (xsec_->is_material_fissile.at(material_id)) {    

    // Fill the two vectors with the appropriate fission sources, summed
    // over all groups
    for (int group_in = 0; group_in < n_group_; ++group_in) {
      // Retrieve cell scalar flux, inverse sigma t, fission terms
      std::vector<double> group_cell_scalar_flux(n_q_);
      this->GetGroupCellScalarFlux(group_cell_scalar_flux, group_in);
      auto inv_sigma_t = xsec_->inv_sigt.at(material_id)[g];
      auto scaled_fission_transfer =
          scaled_fiss_transfer_.at(material_id)(group_in, g);
      
      for (int q = 0; q < n_q_; ++q) {
        cell_q[q] +=
            scaled_fission_transfer * group_cell_scalar_flux[q];
        cell_q_over_total[q] +=
            cell_q[q] * inv_sigma_t;
      }
    }    
  }

  //Integrate and add both source terms
  for (int q = 0; q < n_q_; ++q) {
    
    cell_q[q] *= fv_->JxW(q);
    cell_q_over_total[q] *= fv_->JxW(q);
    
    for (int i = 0; i < dofs_per_cell_; ++i) {
      // First scattering term
      cell_rhs(i) += fv_->shape_value(i, q) * cell_q[q];
      // Second scattering term
      cell_rhs(i) +=
          omega_[dir] * fv_->shape_grad(i, q) * cell_q_over_total[q];
    }
  }
}

template<int dim>
void SelfAdjointAngularFlux<dim>::IntegrateScatteringLinearForm (
      typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      dealii::Vector<double> &cell_rhs,
      const int &g,
      const int &dir) {
  // Get material id for the given cell:
  int material_id = cell->material_id();

  // SAAF has two scattering terms, the first is proportional to scalar flux
  // time sigma_s over 4pi, the second has an additional division by sigma_t.
  // These two vectors hold those values, respectively,  at each quadrature
  // point:
  std::vector<double> cell_scatter_flux(n_q_);
  std::vector<double> cell_scatter_over_total_flux(n_q_);
  
  // Iterate over groups to populate cell_scatter_flux
  for (int group_in = 0; group_in < n_group_; ++group_in) {
    std::vector<double> group_cell_scalar_flux(n_q_);
    this->GetGroupCellScalarFlux(group_cell_scalar_flux, group_in);

    // Get needed cross-sections:
    auto sigma_s_per_ster = xsec_->sigs_per_ster.at(material_id)(group_in, g);
    auto inv_sigma_t = xsec_->inv_sigt.at(material_id)[g];
    
    // Fold group cell scalar flux into total cell scalar flux and multiply by
    // appropriate cross-sections:
    for (int q = 0; q < n_q_; ++q) {
      cell_scatter_flux[q] += sigma_s_per_ster * group_cell_scalar_flux[q];
      cell_scatter_over_total_flux[q] = cell_scatter_flux[q] * inv_sigma_t;
    }
  }

  // Integrate and add both scattering terms
  for (int q = 0; q < n_q_; ++q) {
    cell_scatter_flux[q] *= fv_->JxW(q);
    cell_scatter_over_total_flux[q] *= fv_->JxW(q);
    for (int i = 0; i < dofs_per_cell_; ++i) {
      // First scattering term
      cell_rhs(i) += fv_->shape_value(i, q) * cell_scatter_flux[q];
      // Second scattering term
      cell_rhs(i) +=
          omega_[dir] * fv_->shape_grad(i, q) * cell_scatter_over_total_flux[q];
    }
  }
}

template <int dim>
void SelfAdjointAngularFlux<dim>::PreassembleCellMatrices () {
  // Reinitialize FEM values to an arbitrary cell (the first one) in the list of
  // local cells.
  fv_->reinit(dat_ptr_->local_cells[0]);

  // For each quadrature angle, generate the Collision and Streaming matrices
  for (int q = 0; q < n_q_; ++q) {

    dealii::FullMatrix<double> temp_matrix(dofs_per_cell_, dofs_per_cell_);
    
    for (int i = 0; i < dofs_per_cell_; ++i) {
      for (int j = 0; j < dofs_per_cell_; ++j) {
        temp_matrix(i,j) =
            (fv_->shape_value(i,q) * fv_->shape_value(j,q));
      }
    }
    
    pre_collision_[q] = temp_matrix;

    temp_matrix = 0;
    
    // Streaming terms also depend on direction
    for (int dir = 0; dir < n_dir_; ++dir) {
      for (int i = 0; i < dofs_per_cell_; ++i) {
        for (int j = 0; j < dofs_per_cell_; ++j) {
          temp_matrix(i,j) =
              (fv_->shape_grad(i,q) * omega_[dir])
              *
              (fv_->shape_grad(j,q) * omega_[dir]);
        }
      }
      pre_streaming_[{dir, q}] = temp_matrix;
    }
  }
}

/*
 * =============================================================================
 * PROTECTED FUNCTIONS
 * =============================================================================
 */  

template <int dim>
void SelfAdjointAngularFlux<dim>::GetGroupCellScalarFlux
(std::vector<double> &to_fill, int group) {

  // Get the global scalar flux for the current group
  auto & group_global_scalar_flux =
      mat_vec_->moments[equ_name_][std::make_tuple(group, 0, 0)];
  // Evaluate the global scalar flux at the quadrature points of the current
  // cell and store in return_vector (dealii function in FEValuesBase)
  fv_->get_function_values(group_global_scalar_flux, to_fill);
  
} 

template class SelfAdjointAngularFlux<1>;
template class SelfAdjointAngularFlux<2>;
template class SelfAdjointAngularFlux<3>;

