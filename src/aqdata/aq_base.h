#ifndef __aq_base_h__
#define __aq_base_h__

#include <deal.II/base/tensor.h>
#include <deal.II/base/parameter_handler.h>

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

using namespace dealii;

//! This class provides AQ data, component indices and reflective directions.
/*!
 This class is the base class of AQ data. The directions are represented by a
 series of <a href="https://www.dealii.org/8.5.0/doxygen/deal.II/classTensor.htm
 l" style="color:blue"><b>Tensor<1,dim></b></a>. deal.II has well defined
 operations such as dot product of tensors.
 
 \author Weixiong Zheng
 \date 2017/04
 */
template <int dim>
class AQBase
{
public:
  /*!
   Class constructor.
   
   \param prm ParameterHandler object.
   */
  AQBase<dim> (ParameterHandler &prm);
  
  //! Virtual destructor.
  virtual ~AQBase ();

  /*!
   A virtual function to produce angular quadrature.
   
   \note One has to override this function in derived classes.
   \return Void.
   */
  virtual void produce_angular_quad ();
  virtual void initialize_component_index ();
  
  /*!
   This function calls other members to:
   
   (1) Create aq data;
   
   (2) Create indices of components in transport system;
   
   (3) Create indices of reflective directions for all directions in the aq on all
       boundaries.
   
   \return Void.
   */
  void make_aq ();
  
  /*!
   This function output angular quadrature in aq.txt file.
   
   \return Void.
   */
  void print_angular_quad ();

  /*!
   A function to return SN order in integer.
   
   \return SN order.
   */
  unsigned int get_sn_order ();
  
  /*!
   A function to return total number of directions.
   
   \return Total number of directions in integer.
   */
  unsigned int get_n_dir ();
  
  /*!
   A function to return total number of components in HO equation.
   
   \return Total number of components in HO equation.
   */
  unsigned int get_n_total_ho_vars ();
  
  /*!
   A function to return AQBase<dim>::wi.
   
   \return A vector of all angular weights.
   */
  std::vector<double> get_angular_weights ();
  
  /*!
   A function to return all the directions.
   
   \return A vector of dealii::Tensor<1, dim> representing directions.
   */
  std::vector<Tensor<1, dim> > get_all_directions ();
  
  /*!
   A function to return HO component indices, AQBase<dim>::component_index.
   
   \return A std::map for (group_idx, dir_idx)->component_idx.
   */
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
  get_component_index_map ();
  
  /*!
   A function to return AQBase<dim>::inverse_component_index.
   
   \return A Hash table for component_idx->(group_idx, dir_idx).
   */
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> >
  get_inv_component_map ();
  
  /*!
   A function to return AQBase<dim>::reflective_direction_index.
   
   \return std::map for the mapping: (boundary_id, current_dir)->refl_dir.
   */
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
  get_reflective_direction_index_map ();

protected:
  const double pi;//!< Constant PI, 3.14159...
  double total_angle;//!< Total solid angle. 
  bool have_reflective_bc;//!< Boolean about if there is any reflective boundary.
  std::string transport_model_name;//!< Transport model name in string.
  std::string discretization;//!< Spatial discretization method in string.
  unsigned int n_azi;//!< Total number of azimuthal angles.
  unsigned int n_group;//!< Total number of groups.
  unsigned int n_dir;//!< Total number of directions in the quadrature.
  unsigned int n_total_ho_vars;//!< Total number of components in HO equation.
  std::vector<Tensor<1, dim> > omega_i;//!< All directions in Tensor<1, dim>
  std::vector<double> wi;//!< All angular weights
  
  /*!
   A std::map using pair of group and direction indices as key and component index
   as value.
   */
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
  component_index;
  
  /*!
   A Hash table using component index as key and pair of group and direction indices
   as value.
   */
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> >
  inverse_component_index;
  
  /*!
   A std::map using pair of boundary id and current direction index as key and 
   reflective direction index as value.
   */
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
  reflective_direction_index;

private:
  /*!
   This function returns full name of the angular quadrature in string based on
   the abbreviated name.
   
   \return Full AQ name in string.
   */
  std::string produce_quadrature_name ();
  
  /*!
   This function initialize reflective direction indices on all boundaries. If
   AQBase<dim>::have_reflective_bc is true, this function will generate such indices
   for all boundaries and result in std::map using std::pair of boundary id and 
   current direction index as the key and reflective direction index as the value.
   
   \return Void.
   */
  void initialize_ref_bc_index ();

  std::string aq_name;//! Abbreviated lower-case name of AQ in string.
};

#endif //__aq_base_h__
