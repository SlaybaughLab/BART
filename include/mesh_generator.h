#ifndef __mesh_generator_h__
#define __mesh_generator_h__

#include "problem_definition.h"

template <int dim>
class MeshGenerator
{
public:
  MeshGenerator (const std_cxx11::shared_ptr<ProblemDefinition<dim> > p_def);
  ~MeshGenerator ();
  
  void make_grid (Triangulation<dim> &tria);
  void get_relevant_cell_iterators
  (DoFHandler<dim> &dof_handler,
   std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells,
   std::vector<bool> &is_cell_at_bd,
   std::vector<bool> &is_cell_at_ref_bd);
  
private:
  void generate_initial_grid (Triangulation<dim> &tria);
  void initialize_material_id (Triangulation<dim> &tria);
  void setup_boundary_ids (Triangulation<dim> &tria);
  
  // utility member functions
  void get_cell_relative_position
  (Point<dim> &position,
   std::vector<unsigned int> &relative_position);
  
  
  bool is_mesh_generated;
  bool have_reflective_bc;
  std::string mesh_filename;
  unsigned int global_refinements;
  std::map<std::vector<unsigned int>, unsigned int> relative_position_to_id;
  std::unordered_map<unsigned int, bool> is_reflective_bc;
  
  std::vector<double> axis_max_values;
  std::vector<double> cell_size_all_dir;
  std::vector<unsigned int> ncell_per_dir;
};

#endif //__mesh_generator_h__
