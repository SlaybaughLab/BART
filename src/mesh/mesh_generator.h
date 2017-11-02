#ifndef __mesh_generator_h__
#define __mesh_generator_h__

#include <deal.II/distributed/tria.h>
#include <deal.II/base/parameter_handler.h>

#include <unordered_map>
#include <map>
#include <vector>

using namespace dealii;

//! This class provides functionalities to generate a distributed mesh.
/*!
 This class implement generating meshes using user-defined parameters. Supported
 functionalities are:
 
 (1) Genereate a coarse mesh;
 
 (2) Set up material IDs for all cells in coarse mesh;
 
 (3) Perform global refinements to the mesh.
 
 \author Weixiong Zheng
 \date 2017/05
 */
template <int dim>
class MeshGenerator
{
public:
  MeshGenerator (ParameterHandler &prm);
  ~MeshGenerator ();
  
  /*!
   This function generates or reads in coarse mesh without global refinements.
   Currently, read-in is not fully implemented yet. 
   
   \note There is unknown error if total number of cells in the coarse mesh 
   cannot be divided by number of processors.
   
   \todo Add functionality to read in meshes.
   
   \param tria Triangulation object.
   \return Void. Modify tria in place.
   */
  void make_grid (parallel::distributed::Triangulation<dim> &tria);
  
  /*!
   This function initializes iterators for cells on current processor.
   
   \param dof_handler An DoFHandler object containing iterators of all cells.
   \param local_cells Vector of active cell iterators living only on current 
   processor.
   \return Void.
   */
  void get_relevant_cell_iterators
  (const DoFHandler<dim> &dof_handler,
   std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells);
  
  /*!
   Public member function to get total number of global refinements.
   
   \return Total number of global refinements.
   */
  unsigned int get_uniform_refinement ();
  
  /*!
   Public member function to get Hash table showing if a boundary is reflective.
   
   \return std::unordered_map with key as the boundary_id (integer) and value as
   the boolean.
   */
  std::unordered_map<unsigned int, bool> get_reflective_bc_map ();
  
private:
  
  /*!
   Generate initial coarse grid according to user defined parameters. The mesh
   is a hyer rectangle. In 2D, it is a rectangle. In 3D it is a cuboid. Defining
   the hyper rectangle requires two diagonal points and mesh cells per axis.
   
   \param tria Triangulation object.
   \return Void. Modify tria in place.
   */
  void generate_initial_grid (parallel::distributed::Triangulation<dim> &tria);
  
  /*!
   This member function set up material IDs to the cells belonging to current
   processor on the coarse mesh before performing global refinements.
   
   \param tria Triangulation object.
   \return Void. Modify tria in place.
   */
  void initialize_material_id (parallel::distributed::Triangulation<dim> &tria);
  
  /*!
   This function set up boundary IDs. The naming philosophy is xmin->0, xmax->1,
   ymin->2, ymax->3, zmin->4, zmax->5 if boundaries are applicable.
   
   \param tria Triangulation object.
   \return Void. Modify tria in place.
   */
  void setup_boundary_ids (parallel::distributed::Triangulation<dim> &tria);
  void initialize_relative_position_to_id_map (ParameterHandler &prm);
  void preprocess_reflective_bc (ParameterHandler &prm);
  void process_coordinate_information (ParameterHandler &prm);
  // utility member functions
  void get_cell_relative_position
  (Point<dim> &position,
   std::vector<unsigned int> &relative_position);
  
  /*!
   Boolean to determine if mesh needs to be generated or read-in. Currently, BART
   can only use generated mesh and read-in functionality is not fully developed.
   */
  bool is_mesh_generated;
  bool have_reflective_bc;//!< Boolean to determine if reflective BC is used.
  std::string mesh_filename;//!< Mesh filename if mesh is read in.
  unsigned int global_refinements;//!< Number of global refinements to perform.
  std::map<std::vector<unsigned int>, unsigned int> relative_position_to_id;
  std::unordered_map<unsigned int, bool> is_reflective_bc;
  
  std::vector<double> axis_max_values;//!< Max values per axis in the mesh.
  std::vector<double> cell_size_all_dir;//!< Cell length per direction on the coarse mesh.
  std::vector<unsigned int> ncell_per_dir;//!< Initial number of cells per axis on the coarse mesh.
};

#endif //__mesh_generator_h__
