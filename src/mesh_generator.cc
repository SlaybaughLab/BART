#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>

#include <fstream>

#include "../include/mesh_generator.h"

template <int dim>
MeshGenerator<dim>::MeshGenerator
(std_cxx11::shared_ptr<ProblemDefinition<dim> > p_def)
{
  axis_max_values = p_def->get_axis_maxes ();
  ncell_per_dir = p_def->get_ncells ();
  cell_size_all_dir = p_def->get_cell_sizes ();
  is_mesh_generated = p_def->get_generated_mesh_bool ();
  if (!is_mesh_generated)
    mesh_filename = p_def->get_mesh_file_name ();
  global_refinements = p_def->get_uniform_refinement ();
  relative_position_to_id = p_def->get_id_map ();
  have_reflective_bc = p_def->get_reflective_bool ();
  if (have_reflective_bc)
    is_reflective_bc = p_def->get_reflective_bc_map ();
}

template <int dim>
MeshGenerator<dim>::~MeshGenerator ()
{
}

template <int dim>
void MeshGenerator<dim>::make_grid (Triangulation<dim> &tria)
{
  std::cout << "generated mesh? " << is_mesh_generated << std::endl;
  if (is_mesh_generated)
  {
    generate_initial_grid (tria);
    initialize_material_id (tria);
    setup_boundary_ids (tria);
    tria.refine_global (global_refinements);
  }
  else
  {
    GridIn<dim> gi;
    gi.attach_triangulation (tria);
    std::ifstream f(mesh_filename);
    gi.read_msh (f);
  }
}

template <int dim>
void MeshGenerator<dim>::generate_initial_grid (Triangulation<dim> &tria)
{
  std::cout << "make initial grid" << std::endl;
  Point<dim> origin;
  Point<dim> diagonal;
  switch (dim)
  {
    case 1:
    {
      diagonal[0] = axis_max_values[0];
      break;
    }
      
    case 2:
    {
      diagonal[0] = axis_max_values[0];
      diagonal[1] = axis_max_values[1];
      break;
    }
      
    case 3:
    {
      diagonal[0] = axis_max_values[0];
      diagonal[1] = axis_max_values[1];
      diagonal[2] = axis_max_values[2];
      break;
    }
      
    default:
      break;
  }
  GridGenerator::subdivided_hyper_rectangle (tria,
                                             ncell_per_dir,
                                             origin,
                                             diagonal);
}

template <int dim>
void MeshGenerator<dim>::initialize_material_id (Triangulation<dim> &tria)
{
  std::cout << "set up material id" << std::endl;
  AssertThrow (is_mesh_generated==true,
               ExcMessage("mesh read in have to have boundary ids associated"));
  for (auto cell=tria.begin_active();
       cell!=tria.end(); ++cell)
    if (cell->is_locally_owned())
    {
      Point<dim> center = cell->center ();
      std::vector<unsigned int> relative_position (3);
      get_cell_relative_position (center, relative_position);
      unsigned int material_id = relative_position_to_id[relative_position];
      cell->set_material_id (material_id);
    }
}

template <int dim>
void MeshGenerator<dim>::get_cell_relative_position (Point<dim> &center,
                                                     std::vector<unsigned int> &relative_position)
{
  AssertThrow (relative_position.size()==3,
               ExcMessage("relative position should be size 3 for any dimension"));
  if (dim>=1)
  {
    relative_position[0] = static_cast<unsigned int>(center[0] / cell_size_all_dir[0]);
    if (dim>=2)
    {
      relative_position[1] = static_cast<unsigned int>(center[1] / cell_size_all_dir[1]);
      if (dim==3)
        relative_position[2] = static_cast<unsigned int>(center[2] / cell_size_all_dir[2]);
    }
  }
}

template <int dim>
void MeshGenerator<dim>::setup_boundary_ids (Triangulation<dim> &tria)
{
  std::cout << "set up boundary id" << std::endl;
  AssertThrow (is_mesh_generated==true,
               ExcMessage("mesh read in have to have boundary ids associated"));
  AssertThrow (axis_max_values.size()==dim,
               ExcMessage("number of entries axis max values should be dimension"));
  
  for (auto cell=tria.begin_active();
       cell!=tria.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
      {
        if (cell->face(fn)->at_boundary())
        {
          Point<dim> ct = cell->face(fn)->center();
          // left boundary
          if (std::fabs(ct[0])<1.0e-14)
            cell->face(fn)->set_boundary_id (0);
          
          // right boundary
          if (std::fabs(ct[0]-axis_max_values[0])<1.0e-14)
            cell->face(fn)->set_boundary_id (1);
          
          // 2D and 3D boundaries
          if (dim>1)
          {
            // 2D boundaries
            // front boundary
            if (std::fabs(ct[1])<1.0e-14)
              cell->face(fn)->set_boundary_id (2);
            
            // rear boundary
            if (std::fabs(ct[1]-axis_max_values[1])<1.0e-14)
              cell->face(fn)->set_boundary_id (3);
            
            // 3D boundaries
            if (dim>2)
            {
              // front boundary
              if (std::fabs(ct[2])<1.0e-14)
                cell->face(fn)->set_boundary_id (4);
              
              // rear boundary
              if (std::fabs(ct[2]-axis_max_values[2])<1.0e-14)
                cell->face(fn)->set_boundary_id (5);
            }
          }
        }
      }// face
    }// locally owned cell
  }// cell
}

template <int dim>
void MeshGenerator<dim>::get_relevant_cell_iterators
(DoFHandler<dim> &dof_handler,
 std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells,
 std::vector<bool> &is_cell_at_bd,
 std::vector<bool> &is_cell_at_ref_bd)
{
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    if (cell->is_locally_owned())
    {
      local_cells.push_back (cell);
      if (cell->at_boundary())
        is_cell_at_bd.push_back (true);
      else
        is_cell_at_bd.push_back (false);
      
      if (have_reflective_bc)
      {
        bool at_ref_bd = false;
        for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
          if (cell->at_boundary(fn))
          {
            unsigned int bd_id = cell->face(fn)->boundary_id ();
            if (!at_ref_bd && is_reflective_bc[bd_id])
              at_ref_bd = true;
          }
        if (at_ref_bd)
          is_cell_at_ref_bd.push_back (true);
        else
          is_cell_at_ref_bd.push_back (false);
      }
    }
}

template class MeshGenerator<2>;
template class MeshGenerator<3>;
