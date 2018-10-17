#include "self_adjoint_angular_flux.h"

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
    std::fill(cell_q.begin(), cell_q.end(), q_per_ster);
    std::fill(cell_q_over_total.begin(), cell_q_over_total.end(),
              q_per_ster * xsec_->inv_sigt.at(material_id)[g]);
    
  } else if (xsec_->is_material_fissile.at(material_id)) {    
    
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
    std::vector<double> group_cell_scalar_flux =
        this->GetGroupScalarFlux(group_in);

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
    pre_collision_[q] = this->CellCollisionMatrix(q);
    // Streaming terms also depend on direction
    for (int dir = 0; dir < n_dir_; ++dir)
      pre_streaming_[{dir, q}] = this->CellStreamingMatrix(q, dir);
  }
}

/*
 * =============================================================================
 * PROTECTED FUNCTIONS
 * =============================================================================
 */  

template <int dim>
dealii::FullMatrix<double>
SelfAdjointAngularFlux<dim>::CellCollisionMatrix (int q) {
  
  dealii::FullMatrix<double> return_matrix (dofs_per_cell_, dofs_per_cell_);
  
  for (int i = 0; i < dofs_per_cell_; ++i) {
      for (int j = 0; j < dofs_per_cell_; ++j) {
        return_matrix(i,j) =
            (fv_->shape_value(i,q) * fv_->shape_value(j,q));
      }
  }
  
  return return_matrix;
}

template <int dim>
dealii::FullMatrix<double>
SelfAdjointAngularFlux<dim>::CellStreamingMatrix (int q, int dir) {
  
  dealii::FullMatrix<double> return_matrix (dofs_per_cell_, dofs_per_cell_);
  
  for (int i = 0; i < dofs_per_cell_; ++i) {
      for (int j = 0; j < dofs_per_cell_; ++j) {
        return_matrix(i,j) =
            (fv_->shape_grad(i,q) * omega_[dir])
            *
            (fv_->shape_grad(j,q) * omega_[dir]);
      }
  }
  
  return return_matrix;
}

template <int dim>
std::vector<double>
SelfAdjointAngularFlux<dim>::GetGroupCellScalarFlux(int group) {
  std::vector<double> return_vector(n_q_);

  // Get the global scalar flux for the current group
  auto group_global_scalar_flux =
      mat_vec_->moments[equ_name_][std::make_tuple(group, 0, 0)];
  // Evaluate the global scalar flux at the quadrature points of the current
  // cell and store in return_vector (dealii function in FEValuesBase)
  fv_->get_function_values(group_global_scalar_flux, return_vector);
  
  return return_vector;
}



