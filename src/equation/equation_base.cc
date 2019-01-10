#include "equation_base.h"

template <int dim>
EquationBase<dim>::EquationBase (
    const std::string &equation_name,
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    equ_name_(equation_name),
    dat_ptr_(dat_ptr),
    mat_vec_(dat_ptr_->mat_vec),
    xsec_(dat_ptr_->xsec),
    is_eigen_problem_(prm.get_bool("do eigenvalue calculations")),
    do_nda_(prm.get_bool("do nda")),
    have_reflective_bc_(prm.get_bool("have reflective boundary")),
    p_order_(dat_ptr_->fe_data.p_order[equ_name_]),
    n_dir_((equ_name_=="nda"||equ_name_=="tg_nda"||equ_name_=="diffusion") ? 1 :
        dat_ptr_->aq->GetNDir()),
    n_group_(prm.get_integer("number of groups")),
    n_total_vars_(equ_name_=="tg_nda"?1:n_group_*n_dir_),
    discretization_(dat_ptr_->fe_data.discretization[equ_name_]),
    is_reflective_bc_(dat_ptr_->mesh.GetReflectiveBCMap()),
    fv_(dat_ptr_->fe_data.fv[equ_name_]),
    fvf_(dat_ptr_->fe_data.fvf[equ_name_]),
    fvf_nei_(dat_ptr_->fe_data.fvf_nei[equ_name_]),
    n_q_(dat_ptr_->fe_data.n_q[equ_name_]),
    n_qf_(dat_ptr_->fe_data.n_qf[equ_name_]),
    dofs_per_cell_(dat_ptr_->fe_data.dofs_per_cell[equ_name_]),
    local_dof_indices_(
        std::vector<dealii::types::global_dof_index>(dofs_per_cell_)),
    neigh_dof_indices_(
        std::vector<dealii::types::global_dof_index>(dofs_per_cell_)),
    lin_alg_(prm, equ_name_, n_total_vars_) {
  ProcessInput();
}

template <int dim>
EquationBase<dim>::~EquationBase () {}

template <int dim>
void EquationBase<dim>::ProcessInput () {
  // aq data related
  if (equ_name_!="tg_nda" && equ_name_!="diffusion") {
    ho_comp_ind_ = dat_ptr_->aq->GetCompInd();
    ho_inv_comp_ind_ = dat_ptr_->aq->GetInvCompInd();
    w_ = dat_ptr_->aq->GetAQWeights ();
    omega_ = dat_ptr_->aq->GetAQDirs ();

    if (have_reflective_bc_ && equ_name_!="nda")
      ref_dir_ind_ = dat_ptr_->aq->GetRefDirInd ();
  }
}

// TODO: If NDA is developed, modify this function to contain functionality for
// closure terms
template <int dim>
void EquationBase<dim>::AssembleBilinearForms() {
  dat_ptr_->pcout << "Assemble bilinear forms for " << equ_name_ << std::endl;

  dat_ptr_->pcout << "Assemble volumetric bilinear forms" << std::endl;
  AssembleVolumeBoundaryBilinearForms();

  // interface terms for dfem, fv and rtk
  if (discretization_!="cfem") {
    dat_ptr_->pcout << "Assemble interface bilinear forms for " << equ_name_
        << std::endl;
    AssembleInterfaceBilinearForms();
  }

  // initialize preconditioners for current equation system
  lin_alg_.InitPrecond(
      mat_vec_->sys_mats.at(equ_name_),
      mat_vec_->sys_rhses.at(equ_name_));
}

template <int dim>
void EquationBase<dim>::AssembleVolumeBoundaryBilinearForms () {
  // do the preassembling
  PreassembleCellMatrices();
  dealii::FullMatrix<double> local_mat(dofs_per_cell_, dofs_per_cell_);

  for (int k=0; k<n_total_vars_; ++k) {
    *(mat_vec_->sys_mats[equ_name_][k]) = 0.0;
    int g=GetCompGrpInd(k);
    int dir=GetCompDirInd(k);
    dat_ptr_->pcout
        << "Assemble volume and boundary bilinear forms for Direction " << dir
        << " in Group " << g << std::endl;

    for (auto& cell : dat_ptr_->local_cells) {
      fv_->reinit(cell);
      cell->get_dof_indices(local_dof_indices_);
      local_mat = 0.0;
      // integrate cell bilinear form
      IntegrateCellBilinearForm(cell, local_mat, g, dir);
      for (unsigned int fn=0; fn<dealii::GeometryInfo<dim>::faces_per_cell; ++fn)
        // if cell has face at boundary, do integration
        if (cell->at_boundary(fn)) {
          fvf_->reinit(cell, fn);
          IntegrateBoundaryBilinearForm(cell, fn, local_mat, g, dir);
        }
      // add the local matrix to global matrix
      mat_vec_->sys_mats[equ_name_][k]->add(
          local_dof_indices_,
          local_dof_indices_,
          local_mat);
    }
    mat_vec_
        ->sys_mats[equ_name_][k]->compress(dealii::VectorOperation::add);
  }
}

template <int dim>
void EquationBase<dim>::AssembleInterfaceBilinearForms () {
  dealii::FullMatrix<double> vi_ui (dofs_per_cell_, dofs_per_cell_);
  dealii::FullMatrix<double> vi_ue (dofs_per_cell_, dofs_per_cell_);
  dealii::FullMatrix<double> ve_ui (dofs_per_cell_, dofs_per_cell_);
  dealii::FullMatrix<double> ve_ue (dofs_per_cell_, dofs_per_cell_);

  for (int k=0; k<n_total_vars_; ++k) {
    int g = GetCompGrpInd (k);
    int dir = GetCompDirInd (k);
    if (equ_name_=="nda" || equ_name_=="diffusion")
      dat_ptr_->pcout << "Assemble interface bilinear form for " << equ_name_
          << " in Group " << g << std::endl;
    else
      dat_ptr_->pcout << "Assemble interface bilinear form for " << equ_name_
          << " in Group " << g << ", Direction " << dir << std::endl;
    for (auto& cell : dat_ptr_->local_cells) {
      cell->get_dof_indices (local_dof_indices_);
      for (unsigned int fn=0; fn<dealii::GeometryInfo<dim>::faces_per_cell; ++fn)
        if (!cell->at_boundary(fn) &&
            cell->neighbor(fn)->id()<cell->id()) {
          fvf_->reinit (cell, fn);
          typename dealii::DoFHandler<dim>::cell_iterator
          neigh = cell->neighbor(fn);
          neigh->get_dof_indices (neigh_dof_indices_);
          // note that the same face on the neighboring cell has different index
          // number. So instead, we retrieve the index of current face in
          // neighboring cell using cell->neighbor_face_no(fn)
          fvf_nei_->reinit (neigh, cell->neighbor_face_no(fn));

          vi_ui = 0;
          vi_ue = 0;
          ve_ui = 0;
          ve_ue = 0;

          IntegrateInterfaceBilinearForm (cell, neigh, fn,
                                          vi_ui, vi_ue, ve_ui, ve_ue,
                                          g, dir);
          mat_vec_->sys_mats[equ_name_][k]->add (
              local_dof_indices_,
              local_dof_indices_,
              vi_ui);

          mat_vec_->sys_mats[equ_name_][k]->add (
              local_dof_indices_,
              neigh_dof_indices_,
              vi_ue);

          mat_vec_->sys_mats[equ_name_][k]->add (
              neigh_dof_indices_,
              local_dof_indices_,
              ve_ui);

          mat_vec_->sys_mats[equ_name_][k]->add (
              neigh_dof_indices_,
              neigh_dof_indices_,
              ve_ue);
        }// target faces
    }
    mat_vec_->sys_mats[equ_name_][k]->compress(dealii::VectorOperation::add);
  }// component
}

template <int dim>
void EquationBase<dim>::PreassembleCellMatrices () {}

template <int dim>
void EquationBase<dim>::AssembleFixedLinearForms() {
  for (int k=0; k<n_total_vars_; ++k) {
    int g = GetCompGrpInd (k);
    int dir = GetCompDirInd (k);
    *(mat_vec_->sys_fixed_rhses[equ_name_][k]) = 0.0;
    dealii::Vector<double> cell_rhs (dofs_per_cell_);
    
    for (auto& cell : dat_ptr_->local_cells) {
      cell_rhs = 0.0;
      cell->get_dof_indices (local_dof_indices_);
      fv_->reinit (cell);
      // do cellwise integration
      IntegrateCellFixedLinearForm (cell, cell_rhs, g, dir);
      // map local assembled vector to global
      mat_vec_
          ->sys_fixed_rhses[equ_name_][k]->add(local_dof_indices_, cell_rhs);
    }
    mat_vec_
        ->sys_fixed_rhses[equ_name_][k]->compress (dealii::VectorOperation::add);
  }
}

template <int dim>
void EquationBase<dim>::AssembleLinearForms (const int &g) {
  dealii::Vector<double> cell_rhs (dofs_per_cell_);
  for (int k=0; k<n_total_vars_; ++k)
    if (GetCompGrpInd(k)==g) {
      int dir = GetCompDirInd (k);
      // initialize rhs with fixed source
      *(mat_vec_->sys_rhses[equ_name_][k]) =
          *(mat_vec_->sys_fixed_rhses[equ_name_][k]);

      // cellwise integration
      for (auto& cell : dat_ptr_->local_cells) {
        cell_rhs = 0.0;
        cell->get_dof_indices (local_dof_indices_);
        fv_->reinit (cell);
        // integrate scattering linear form
        IntegrateScatteringLinearForm (cell, cell_rhs, g, dir);
        for (unsigned int fn=0;
             fn<dealii::GeometryInfo<dim>::faces_per_cell; ++fn)
          // integrate boundary linear form at boundary faces
          if (cell->at_boundary(fn)) {
            fvf_->reinit (cell, fn);
            IntegrateBoundaryLinearForm (cell, fn, cell_rhs, g, dir);
          }
        mat_vec_
            ->sys_rhses[equ_name_][k]->add (local_dof_indices_, cell_rhs);
      }
      mat_vec_
          ->sys_rhses[equ_name_][k]->compress (dealii::VectorOperation::add);
    }
}

template <int dim>
void EquationBase<dim>::SolveInGroup (const int &g) {
  // loop over all the components and check corresponding group numbers. Once
  // found, call the linear solvers to solve the equations.

  // Note: redesign is needed in case Krylov method;
  // Overriding could be used when PN-like system is involved
  for (int i=0; i<n_total_vars_; ++i)
    if (GetCompGrpInd(i)==g)
      lin_alg_.LinAlgSolve (
          mat_vec_->sys_mats[equ_name_],
          mat_vec_->sys_flxes[equ_name_],
          mat_vec_->sys_rhses[equ_name_],
          mat_vec_->constraints[equ_name_],
          i);
}

template <int dim>
void EquationBase<dim>::GenerateMoments (
    std::map<std::tuple<int,int,int>, dealii::Vector<double>> &moments,
    std::map<std::tuple<int,int,int>, dealii::Vector<double>> &moments_prev,
    const int &g) {
  auto key = std::make_tuple(g,0,0);
  moments_prev[key] = moments[key];
  // generate moments
  moments[key] = 0.0;
  // note that we only implement scalar flux here
  // TODO: support all other moments
  dealii::Vector<double> aflx_this_proc;
  for (int i=0; i<n_dir_; ++i) {
    aflx_this_proc = *mat_vec_->sys_flxes[equ_name_][i];
    moments[key].add (w_[i], aflx_this_proc);
  }
}

template <int dim>
inline void EquationBase<dim>::ScaleFissTransferMatrices(double keff) {
  if (equ_name_=="diffusion" || equ_name_=="nda")
    std::for_each(xsec_->fiss_transfer.begin(),
        xsec_->fiss_transfer.end(),
        [&](auto &p){
          scaled_fiss_transfer_[p.first] = p.second;
          scaled_fiss_transfer_[p.first] /= keff;
        });
  else
    std::for_each(xsec_->fiss_transfer_per_ster.begin(),
        xsec_->fiss_transfer_per_ster.end(),
        [&](auto &p){
          scaled_fiss_transfer_[p.first] = p.second;
          scaled_fiss_transfer_[p.first] /= keff;
        });
}

template <int dim>
double EquationBase<dim>::EstimateFissSrc () {
  double fiss_src = 0.0;
  for (int g=0; g<n_group_; ++g) {
    // first, estimate local fission source
    
    for (auto& cell : dat_ptr_->local_cells) {
      std::vector<std::vector<double>> local_phis (n_group_,
          std::vector<double> (n_q_));
      // do integration in fissile materials
      int material_id = cell->material_id();
      if (xsec_->is_material_fissile.at(material_id)) {
        fv_->reinit (cell);
        double nu_sigf = xsec_->nu_sigf.at(material_id)[g];
        // retrive solution values scalar flux for current group
        fv_->get_function_values (
            mat_vec_->moments[equ_name_][std::make_tuple(g,0,0)],
            local_phis[g]);
        for (int qi=0; qi<n_q_; ++qi)
          fiss_src += (nu_sigf *
                       local_phis[g][qi] *
                       fv_->JxW(qi));
      }
    }
  }
  // then, we need to accumulate fission source from other processors as well
  return dealii::Utilities::MPI::sum (fiss_src, MPI_COMM_WORLD);
}

template <int dim>
int EquationBase<dim>::GetCompInd (const int &g, const int &dir) const {
  return (equ_name_=="nda"||equ_name_=="diffusion") ? g :
      (equ_name_=="tg_nda") ? 1 : ho_comp_ind_.at({g, dir});
}

template <int dim>
int EquationBase<dim>::GetCompDirInd (const int &comp) const {
  return (equ_name_=="nda"||equ_name_=="diffusion"||equ_name_=="tg_nda") ? 0 :
      ho_inv_comp_ind_.at(comp).second;
}

template <int dim>
int EquationBase<dim>::GetCompGrpInd (const int &comp) const {
  return (equ_name_=="nda"||equ_name_=="diffusion") ? comp :
      (equ_name_=="tg_nda" ? 0 : ho_inv_comp_ind_.at(comp).first);
}

template <int dim>
std::string EquationBase<dim>::GetEquName () const {
  return equ_name_;
}

template <int dim>
int EquationBase<dim>::GetNTotalVars () const {
  return n_total_vars_;
}

template class EquationBase<1>;
template class EquationBase<2>;
template class EquationBase<3>;
