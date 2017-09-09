#include <deal.II/fe/fe_values.h>

#include <boost/algorithm/string.hpp>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "equation_base.h"
#include "bart_tools.h"

using namespace dealii;

template <int dim>
EquationBase<dim>::EquationBase
(std::string equation_name,
 const ParameterHandler &prm,
 const DoFHandler<dim> &dof_handler
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
:
equation_name(equation_name),
n_material(prm.get_integer("number of materials")),
discretization(prm.get("spatial discretization")),
n_group(prm.get_integer("number of groups")),
is_eigen_problem(prm.get_bool("do eigenvalue calculations")),
do_nda(prm.get_bool("do NDA")),
have_reflective_bc(prm.get_bool("have reflective BC")),
p_order(prm.get_integer("finite element polynomial degree")),
nda_quadrature_order(p_order+3) //this is hard coded
{
  // process input for mesh, AQ and material related data
  process_input (msh_ptr, aqd_ptr, mat_ptr);
  // get cell iterators and booleans to tell if they are at boundary on this process
  msh_ptr->get_relevant_cell_iterators (dof_handler, local_cells, is_cell_at_bd);
  alg_ptr = build_linalg (prm, equation_name, n_total_vars)
}

template <int dim>
EquationBase<dim>::~EquationBase ()
{
}

template <int dim>
void EquationBase<dim>::process_input
(const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
{
  // mesh related
  {
    relative_position_to_id = msh_ptr->get_id_map ();
    if (have_reflective_bc)
      is_reflective_bc = msh_ptr->get_reflective_bc_map ();
  }
  
  // aq data related
  {
    // note that n_total_vars will be have to be re-init
    // in derived class of if it's for NDA
    if (equation_name!="nda")
      n_total_vars = aqd_ptr->get_n_total_ho_vars ();
    else
      n_total_vars = n_group;
    n_dir = aqd_ptr->get_n_dir ();
    component_index = aqd_ptr->get_component_index_map ();
    inverse_component_index = aqd_ptr->get_inv_component_map ();
    wi = aqd_ptr->get_angular_weights ();
    omega_i = aqd_ptr->get_all_directions ();
    
    if (have_reflective_bc)
      reflective_direction_index = aqd_ptr->get_reflective_direction_index_map ();
  }
  
  // material related
  {
    all_sigt = mat_ptr->get_sigma_t ();
    all_inv_sigt = mat_ptr->get_inv_sigma_t ();
    all_sigs = mat_ptr->get_sigma_s ();
    all_sigs_per_ster = mat_ptr->get_sigma_s_per_ster ();
    if (is_eigen_problem)
    {
      is_material_fissile = mat_ptr->get_fissile_id_map ();
      all_nusigf = mat_ptr->get_nusigf ();
      all_ksi_nusigf = mat_ptr->get_ksi_nusigf ();
      all_ksi_nusigf_per_ster = mat_ptr->get_ksi_nusigf_per_ster ();
    }
    else
    {
      all_q = mat_ptr->get_q ();
      all_q_per_ster = mat_ptr->get_q_per_ster ();
    }
  }
}

template <int dim>
void EquationBase<dim>::initialize_cell_iterators_this_proc
(const DoFHandler<dim> &dof_handler)
{
  msh_ptr->get_relevant_cell_iterators (dof_handler,
                                        local_cells,
                                        is_cell_at_bd);
}

template <int dim>
void EquationBase<dim>::initialize_system_matrices_vectors
(SparsityPatternType &dsp, IndexSet &local_dofs)
{
  for (unsigned int i=0; i<n_total_vars; ++i)
  {
    sys_mats.push_back (new PETScWrappers::MPI::SparseMatrix);
    sys_mats[i]->reinit (local_dofs,
                         local_dofs,
                         dsp,
                         MPI_COMM_WORLD);
    sys_aflxes.push_back (new PETScWrappers::MPI::Vector);
    sys_aflxes[i]->reinit (local_dofs, MPI_COMM_WORLD);
    sys_rhses.push_back (new PETScWrappers::MPI::Vector);
    sys_rhses[i]->reinit (local_dofs, MPI_COMM_WORLD);
    sys_fixed_rhses.push_back (new PETScWrappers::MPI::Vector);
    sys_fixed_rhses[i]->reinit (local_dofs, MPI_COMM_WORLD);
  }
}

template <int dim>
void EquationBase<dim>::initialize_assembly_related_objects
(FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe)
{
  q_rule = std_cxx11::shared_ptr<QGauss<dim> > (new QGauss<dim> (p_order + 1));
  qf_rule = std_cxx11::shared_ptr<QGauss<dim-1> > (new QGauss<dim-1> (p_order + 1));
  
  fv = std_cxx11::shared_ptr<FEValues<dim> >
  (new FEValues<dim> (*fe, *q_rule,
                      update_values | update_gradients |
                      update_quadrature_points |
                      update_JxW_values));

  fvf = std_cxx11::shared_ptr<FEFaceValues<dim> >
  (new FEFaceValues<dim> (*fe, *qf_rule,
                          update_values | update_gradients |
                          update_quadrature_points | update_normal_vectors |
                          update_JxW_values));

  if (discretization=="dfem")
    fvf_nei = std_cxx11::shared_ptr<FEFaceValues<dim> >
    (new FEFaceValues<dim> (*fe, *qf_rule,
                            update_values | update_gradients |
                            update_quadrature_points | update_normal_vectors |
                            update_JxW_values));

  dofs_per_cell = fe->dofs_per_cell;
  n_q = q_rule->size();
  n_qf = qf_rule->size();

  local_dof_indices.resize (dofs_per_cell);
  neigh_dof_indices.resize (dofs_per_cell);

  if (equation_name=="nda")
  {
    qc_rule = std_cxx11::shared_ptr<QGauss<dim> > (new QGauss<dim> (nda_quadrature_order));
    qfc_rule = std_cxx11::shared_ptr<QGauss<dim-1> > (new QGauss<dim-1> (nda_quadrature_order));
    fvc = std_cxx11::shared_ptr<FEValues<dim> >
    (new FEValues<dim> (*fe, *qc_rule,
                        update_values | update_gradients |
                        update_quadrature_points |
                        update_JxW_values));
    
    fvfc = std_cxx11::shared_ptr<FEFaceValues<dim> >
    (new FEFaceValues<dim> (*fe, *qfc_rule,
                            update_values | update_gradients |
                            update_quadrature_points | update_normal_vectors |
                            update_JxW_values));
    n_qc = qc_rule->size ();
    n_qfc = qfc_rule->size ();
  }
}

template <int dim>
void EquationBase<dim>::assemble_system ()
{
  radio ("Assemble volumetric bilinear forms");
  assemble_volume_boundary_bilinear_form ();
  
  if (discretization=="dfem")
  {
    AssertThrow (equation_name=="ep",
                 ExcMessage("DFEM is only implemented for even parity"));
    radio ("Assemble cell interface bilinear forms for DFEM");
    assemble_interface_bilinear_form ();
  }
}

// TODO: work needed for this assembly facility for NDA
template <int dim>
void EquationBase<dim>::assemble_closure_bilinear_form
(std_cxx11::shared_ptr<EquationBase<dim> > ho_equ_ptr)
{
  // the input is pointer to HO equation pointer s.t. we can have estimation of
  // corrections
  AssertThrow (equation_name=="nda",
               ExcMessage("only instance for NDA calls this function"));
}

template <int dim>
void EquationBase<dim>::assemble_volume_boundary_bilinear_form ()
{
  // volumetric pre-assembly matrices
  std::vector<std::vector<FullMatrix<double> > >
  streaming_at_qp (n_q, std::vector<FullMatrix<double> > (n_dir, FullMatrix<double> (dofs_per_cell, dofs_per_cell)));
  
  std::vector<FullMatrix<double> >
  collision_at_qp (n_q, FullMatrix<double>(dofs_per_cell, dofs_per_cell));
  
  // this sector is for pre-assembling streaming and collision matrices at quadrature
  // points
  {
    typename DoFHandler<dim>::active_cell_iterator cell = local_cells[0];
    fv->reinit (cell);
    pre_assemble_cell_matrices (fv, cell, streaming_at_qp, collision_at_qp);
  }
  
  for (unsigned int k=0; k<n_total_vars; ++k)
  {
    unsigned int g = get_component_group (k);
    unsigned int i_dir = get_component_direction (k);
    radio ("Assembling Component",k,"direction",i_dir,"group",g);
    FullMatrix<double> local_mat (dofs_per_cell, dofs_per_cell);
    
    for (unsigned int ic=0; ic<local_cells.size(); ++ic)
    {
      typename DoFHandler<dim>::active_cell_iterator cell = local_cells[ic];
      fv->reinit (cell);
      cell->get_dof_indices (local_dof_indices);
      local_mat = 0;
      integrate_cell_bilinear_form (cell, cell_rhs, g, i_dir);
      
      if (is_cell_at_bd[ic])
        for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
          if (cell->at_boundary(fn))
          {
            fvf->reinit (cell, fn);
            integrate_boundary_bilinear_form (fvf,
                                              cell,
                                              fn,
                                              local_mat,
                                              g, i_dir);
          }
      sys_mats[k]->add (local_dof_indices,
                        local_dof_indices,
                        local_mat);
    }
    sys_mats[k]->compress (VectorOperation::add);
  }// components
}

// The following is a virtual function for integraing cell bilinear form;
// It can be overriden if cell pre-assembly is desirable
template <int dim>
void EquationBase<dim>::
pre_assemble_cell_matrices
(typename DoFHandler<dim>::active_cell_iterator &cell,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp)
{// this is a virtual function
}

// The following is a virtual function for integraing cell bilinear form;
// It must be overriden
template <int dim>
void EquationBase<dim>::integrate_cell_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 FullMatrix<double> &cell_matrix,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp,
 const unsigned int &g,
 const unsigned int &i_dir=0)
{
}

/** \brief Integrator for boundary weak form per boundary face per angular/group
 *
 * The function is a virtual function. For diffusion-like system, i_dir is set
 * to 0 by default.
 */
template <int dim>
void EquationBase<dim>::integrate_boundary_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 unsigned int &fn,/*face number*/
 FullMatrix<double> &cell_matrix,
 const unsigned int &g,
 const unsigned int &i_dir=0)
{// this is a virtual function. Details have to be provided per transport model.
}

/** \brief Right hand side integrator specifically for boundary terms.
 *
 */
template <int dim>
void EquationBase<dim>::integrate_boundary_linear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 unsigned int &fn,/*face number*/
 Vector<double> &cell_rhses,
 const unsigned int &g,
 const unsigned int &i_dir)
{// this is a virtual function. Details might be provided given different models
}

/** \brief Interface weak form assembly driver.
 * Member function used to assemble interface weak forms. The main functionality
 * is to go through all non-boundary interfaces of the cells owned on current
 * processor and assemble the weak form using interface assembler.
 *
 * There is no need to override this function for SN calculations. Yet, for PN,
 * diffusion etc., this function must be overriden to correctly take care of the
 * angular component.
 */
template <int dim>
void EquationBase<dim>::assemble_interface_bilinear_form ()
{
  FullMatrix<double> vi_ui (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> vi_ue (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ve_ui (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> ve_ue (dofs_per_cell, dofs_per_cell);
  
  for (unsigned int k=0; k<n_total_vars; ++k)
  {
    unsigned int g = get_component_group (k);
    unsigned int i_dir = get_component_direction (k);
    
    for (unsigned int ic=0; ic<local_cells.size(); ++ic)
    {
      typename DoFHandler<dim>::active_cell_iterator
      cell = local_cells[ic];
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
        if (!cell->at_boundary(fn) &&
            cell->neighbor(fn)->id()<cell->id())
        {
          fvf->reinit (cell, fn);
          typename DoFHandler<dim>::cell_iterator
          neigh = cell->neighbor(fn);
          neigh->get_dof_indices (neigh_dof_indices);
          fvf_nei->reinit (neigh, cell->neighbor_face_no(fn));
          
          vi_ui = 0;
          vi_ue = 0;
          ve_ui = 0;
          ve_ue = 0;
          
          integrate_interface_bilinear_form (cell, neigh, fn,
                                             vi_ui, vi_ue, ve_ui, ve_ue,
                                             g, i_dir);
          sys_mats[k]->add (local_dof_indices,
                            local_dof_indices,
                            vi_ui);
          
          sys_mats[k]->add (local_dof_indices,
                            neigh_dof_indices,
                            vi_ue);
          
          sys_mats[k]->add (neigh_dof_indices,
                            local_dof_indices,
                            ve_ui);
          
          sys_mats[k]->add (neigh_dof_indices,
                            neigh_dof_indices,
                            ve_ue);
        }// target faces
    }
    sys_mats[k]->compress(VectorOperation::add);
  }// component
}

/** \brief Virtual function for NDA closure assembly.
 * Main functionality is to loop over all local cells and call closure integrators.
 *
 * This needs to be done
 */
template <int dim>
void EquationBase<dim>::assemble_closure_bilinear_form ()
{
}

// function to provide cell correction vectors at cell quadrature points for NDA
// equation
template <int dim>
void EquationBase<dim>::prepare_cell_corrections
(const std::vector<std::vector<Tensor<1, dim> > > &ho_cell_dpsi,
 const std::vector<Tensor<1, dim> > &ho_cell_dphi,
 const std::vector<double> &ho_cell_phi,
 std::vector<Tensor<1, dim> > &cell_corrections)
{
  AssertThrow (equation_name=="nda",
               ExcMessage("cell corrections have to be calculated in NDA equations"));
  // TODO: fill this up
}

// function to provide boundary correction coefficients at face quadrature points
// for NDA equation
template <int dim>
void EquationBase<dim>::prepare_boundary_corrections
(const std::vector<double> &ho_cell_psi
 std::vector<double> &boundary_corrections)
{
  AssertThrow (equation_name=="nda",
               ExcMessage("face corrections have to be calculated in NDA equations"));
  // TODO: fill this up
}

/** \brief Virtual function for interface integrator.
 * When DFEM is used, this function can be overridden as interface weak form
 * assembler per face per angular and group component.
 *
 * When overridden for diffusion calculations, direction component is set to be
 * zero by default
 */
// The following is a virtual function for integrating DG interface for HO system
// it must be overriden
template <int dim>
void EquationBase<dim>::integrate_interface_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 typename DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
 unsigned int &fn,/*concerning face number in local cell*/
 FullMatrix<double> &vi_ui,
 FullMatrix<double> &vi_ue,
 FullMatrix<double> &ve_ui,
 FullMatrix<double> &ve_ue,
 const unsigned int &g,
 const unsigned int &i_dir=0)
{
}

// generate moments on current process for all groups at once
template <int dim>
void EquationBase<dim>::generate_moments
(std::vector<Vector<double> > &sflx_proc,
 std::vector<Vector<double> > &sflx_proc_old)
{
  // TODO: only scalar flux is generated for now, future will be moments considering
  // anisotropic scattering
  for (unsigned int g=0; g<n_group; ++g)
  {
    sflx_proc_old[g] = sflx_proc[g];
    sflx_proc[g] = 0;
    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
      sflx_proc[g].add (wi[i_dir], *sys_aflxes[get_component_index(i_dir, g)]);
  }
}

// generate specific
template <int dim>
void EquationBase<dim>::generate_moments
(Vector<double> &sflx_proc,
 Vector<double> &sflx_proc_old,
 const unsigned int &g)
{
  // TODO: only scalar flux is generated for now, future will be moments considering
  // anisotropic scattering.
  
  // The difference from the other function is that we generate moments for a
  // specific group g
  AssertThrow (equation_name!="nda",
               ExcMessage("NDA is not supposed to call this function"));
  sflx_proc_old = sflx_proc;
  sflx_proc = 0;
  for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    sflx_proc.add (wi[i_dir], *sys_aflxes[get_component_index(i_dir, g)]);
}

// generate scalar flux from HO solver for NDA
template <int dim>
void EquationBase<dim>::generate_ho_sflx
(std::vector<Vector<double> > &ho_aflx_proc,
 std::vector<Vector<double> > &ho_sflx_proc)
{
  AssertThrow (equation_name=="nda",
               ExcMessage("only NDA is supposed to call this function"));
  for (unsigned int g=0; g<n_group; ++g)
  {
    ho_sflx_proc[g].equ (wi[0], ho_aflx_proc[get_component_index(0, g)]);
    for (unsigned int i_dir=1; i_dir<n_dir; ++i_dir)
      ho_sflx_proc[g].add (wi[i_dir], ho_aflx_proc[get_component_index(i_dir, g)]);
  }
}

// get a copy of HO angular flux at local processor
// The intension is s.t. we can use HO angular flux to generate corrections for NDA
template <int dim>
void EquationBase<dim>::get_ho_aflx (std::vector<Vector<double> > &ho_aflx_proc)
{
  AssertThrow (equation_name!="NDA",
               ExcMessage("this function is supposed to be called in HO equation"));
  ho_aflx_proc.resize (n_total_vars);
  for (unsigned int i=0; i<n_total_vars; ++i)
    ho_aflx_proc[i] = *sys_aflxes[i];
}

template <int dim>
void EquationBase<dim>::scale_fiss_transfer_matrices (double keff)
{
  // TODO: after changing ksi_nusigf to std::vector<Vector<double> >, we'll
  // redo this function
  AssertThrow(is_eigen_problem,
              ExcMessage("Only eigen problem calls this member"));
  scaled_fiss_transfer_per_ster.resize (n_material);
  for (unsigned int m=0; m<n_material; ++m)
  {
    std::vector<std::vector<double> > tmp (n_group, std::vector<double>(n_group));
    if (is_material_fissile[m])
      for (unsigned int gin=0; gin<n_group; ++gin)
        for (unsigned int g=0; g<n_group; ++g)
          tmp[gin][g] = all_ksi_nusigf_per_ster[m][gin][g] / keff;
    scaled_fiss_transfer_per_ster[m] = tmp;
  }
}

// generate rhs for equation
template <int dim>
void EquationBase<dim>::assemble_linear_form
(std::vector<Vector<double> > &sflx_this_proc,
 unsigned int &g)
{
  for (unsigned int k=0; k<this->n_total_vars; ++k)
    if (get_component_group(k)==g)
    {
      unsigned int i_dir = get_component_direction (k);
      *sys_rhses[k] = *sys_fixed_rhses[k];
      for (unsigned int ic=0; ic<this->local_cells.size(); ++ic)
      {
        Vector<double> cell_rhs (this->dofs_per_cell);
        typename DoFHandler<dim>::active_cell_iterator cell = local_cells[ic];
        cell->get_dof_indices (this->local_dof_indices);
        fv->reinit (cell);
        std::vector<double> cell_sflx;
        fv->get_function_values (sflx_this_proc[g], cell_sflx);
        integrate_scattering_linear_form (cell, cell_rhs,
                                          g, i_dir);
        if (is_cell_at_bd[ic])
          for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
            if (cell->at_boundary(fn))
            {
              fvf->reinit (cell, fn);
              integrate_boundary_linear_form (cell, fn, cell_rhs,
                                              cell_sflx,
                                              g, i_dir);
            }
        sys_rhses[k]->add (local_dof_indices, cell_rhs);
      }
      sys_rhses[k]->compress (VectorOperation::add);
    }
}

template <int dim>
void EquationBase<dim>::integrate_scattering_linear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 Vector<double> &cell_rhs,
 std::vector<Vector<double> > &sflx_proc,
 const unsigned int &g,
 const unsigned int &i_dir)
{
}

template <int dim>
void EquationBase<dim>::assemble_fixed_linear_form
(std::vector<Vector<double> > &sflx_prev)
{
  for (unsigned int k=0; k<n_total_vars; ++k)
  {
    unsigned int g = get_component_group (k);
    unsigned int i_dir = get_component_direction (k);
    *sys_fixed_rhses[k] = 0.0;
    for (unsigned int ic=0; ic<local_cells.size(); ++ic)
    {
      Vector<double> cell_rhs (this->dofs_per_cell);
      typename DoFHandler<dim>::active_cell_iterator cell = local_cells[ic];
      cell->get_dof_indices (this->local_dof_indices);
      fv->reinit (cell);
      integrate_cell_fixed_linear_form (cell, cell_rhs,
                                        sflx_prev,
                                        g, i_dir);
      sys_fixed_rhses[k]->add (local_dof_indices, cell_rhs);
    }
    sys_fixed_rhses[k]->compress (VectorOperation::add);
  }
}

template <int dim>
void EquationBase<dim>::integrate_cell_fixed_linear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 Vector<double> &cell_rhs,
 std::vector<Vector<double> > &sflx_prev,
 const unsigned int &g,
 const unsigned int &i_dir)
{
}

template <int dim>
void EquationBase<dim>::initialize_preconditioners ()
{
  alg_ptr->initialize_preconditioners (sys_mats, sys_rhses);
}

template <int dim>
void EquationBase<dim>::solve_in_group (const unsigned int &g)
{
  // loop over all the components and check corresponding group numbers. Once
  // found, call the linear solvers to solve the equations.
  
  // Note: redesign is needed in case Krylov method;
  // Overriding could be used when PN-like system is involved
  for (unsigned int i=0; i<n_total_vars; ++i)
    if (get_component_group(i)==g)
      alg_ptr->linear_algebra_solve (sys_mats, sys_aflxes, sys_rhses, i);
}

template <int dim>
double EquationBase<dim>::estimate_fiss_source
(std::vector<Vector<double> > &phis_this_process)
{
  // first, estimate local fission source
  double fiss_source = 0.0;
  for (unsigned int ic=0; ic<local_cells.size(); ++ic)
  {
    typename DoFHandler<dim>::active_cell_iterator cell = local_cells[ic];
    std::vector<std::vector<double> > local_phis (n_group,
                                                  std::vector<double> (n_q));
    unsigned int material_id = cell->material_id ();
    if (is_material_fissile[material_id])
    {
      fv->reinit (cell);
      for (unsigned int g=0; g<n_group; ++g)
        fv->get_function_values (phis_this_process[g],
                                 local_phis[g]);
      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int g=0; g<n_group; ++g)
          fiss_source += (all_nusigf[material_id][g] *
                          local_phis[g][qi] *
                          fv->JxW(qi));
    }
  }
  // then, we need to accumulate fission source from other processors as well
  return Utilities::MPI::sum (fiss_source, MPI_COMM_WORLD);
}

// wrapper functions used to retrieve info from various Hash tables
template <int dim>
unsigned int EquationBase<dim>::get_component_index
(unsigned int incident_angle_index, unsigned int g)
{
  // retrieve component indecis given direction and group
  // must be used after initializing the index map
  return component_index[std::make_pair (incident_angle_index, g)];
}

template <int dim>
unsigned int EquationBase<dim>::get_component_direction (unsigned int comp_ind)
{
  return inverse_component_index[comp_ind].first;
}

template <int dim>
unsigned int EquationBase<dim>::get_component_group (unsigned int comp_ind)
{
  return inverse_component_index[comp_ind].second;
}

template <int dim>
unsigned int EquationBase<dim>::get_reflective_direction_index
(unsigned int boundary_id, unsigned int incident_angle_index)
{
  AssertThrow (is_reflective_bc[boundary_id],
               ExcMessage ("must be reflective boundary to retrieve the reflective boundary"));
  return reflective_direction_index[std::make_pair (boundary_id,
                                                    incident_angle_index)];
}

// explicit instantiation to avoid linking error
template class EquationBase<2>;
template class EquationBase<3>;
