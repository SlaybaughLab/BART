#ifndef BART_SRC_AQDATA_AQ_BASE_H__
#define BART_SRC_AQDATA_AQ_BASE_H__

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include <deal.II/base/tensor.h>
#include <deal.II/base/parameter_handler.h>

//! This class provides angular quadrature related data.
/*!
 This class is the base class of angular quadrature (AQ) data. The directions 
 are represented by a series of <a href="https://www.dealii.org/8.5.0/doxygen/de
 al.II/classTensor.html" style="color:blue"><b>Tensor<1,dim></b></a>. deal.II 
 has well defined operations such as dot product of tensors.
 
 The main functionalities include:
 
 (1) Provide abstract function to produce angular quadrature;
 
 (2) Provide functions to generate component index related data;
 
 (3) Provide functions to generate mapping for retrieving reflective directions.
 
 \note The restarted version of BART uses (group, direction)->component. Before,
 the direction is the first and group is the second of the key pair.
 
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
  AQBase (dealii::ParameterHandler &prm);
  
  //! Virtual destructor.
  virtual ~AQBase ();

  /*!
   A pure virtual function to produce angular quadrature.
   
   Only 1D Gauss-Legendre is implemented in base class. For multi-D, an
   overriding has to be provided per derived class.
   
   \note One has to override this function in derived classes for multi-D or 
   non-Gauss-Legendre 1D angular quadratures.
   \return Void.
   */
  virtual void produce_angular_quad ();
  
  /*!
   A virtual function to initialize component index given group and direction 
   index. By default, this indexing is for SN method.
   */
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
  int get_sn_order ();
  
  /*!
   A function to return total number of directions.
   
   \return Total number of directions in integer.
   */
  int get_n_dir ();
  
  /*!
   A function to return total number of components in HO equation.
   
   \return Total number of components in HO equation.
   */
  int get_n_total_ho_vars ();
  
  /*!
   A function to return AQBase<dim>::wi.
   
   \return A vector of all angular weights.
   */
  std::vector<double> get_angular_weights ();
  
  /*!
   A function to return all the directions.
   
   \return A vector of dealii::Tensor<1, dim> representing directions.
   */
  std::vector<dealii::Tensor<1, dim>> get_all_directions ();
  
  /*!
   A function to return HO component indices, AQBase<dim>::component_index.
   
   \return A std::map for (group_idx, dir_idx)->component_idx.
   */
  std::map<std::pair<int, int>, int> get_component_index_map ();
  
  /*!
   A function to return AQBase<dim>::inverse_component_index_.
   
   \return A Hash table for component_idx->(group_idx, dir_idx).
   */
  std::unordered_map<int, std::pair<int, int>> get_inv_component_map ();
  
  /*!
   A function to return AQBase<dim>::reflective_direction_index.
   
   \return std::map for the mapping: (boundary_id, current_dir)->refl_dir.
   */
  std::map<std::pair<int, int>, int> get_reflective_direction_index_map ();

protected:
  const double k_pi;//!< Constant PI, 3.14159...
  double total_angle_;//!< Total solid angle.
  bool have_reflective_bc_;//!< Boolean about if there is any reflective boundary.
  std::string transport_model_name_;//!< Transport model name in string.
  std::string discretization_;//!< Spatial discretization method in string.
  int n_azi_;//!< Total number of azimuthal angles.
  int n_group_;//!< Total number of groups.
  int n_dir_;//!< Total number of directions in the quadrature.
  int n_total_ho_vars_;//!< Total number of components in HO equation.
  std::vector<dealii::Tensor<1, dim>> omega_i_;//!< All directions in Tensor<1, dim>
  std::vector<double> wi_;//!< All angular weights
  
  /*!
   A std::map using pair of group and direction indices as key and component index
   as value.
   */
  std::map<std::pair<int, int>, int> component_index_;
  
  /*!
   A Hash table using component index as key and pair of group and direction indices
   as value.
   */
  std::unordered_map<int, std::pair<int, int>> inverse_component_index_;
  
  /*!
   A std::map using pair of boundary id and current direction index as key and
   reflective direction index as value.
   */
  std::map<std::pair<int, int>, int> reflective_direction_index_;
  
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
  
  std::string aq_name_;//! Abbreviated lower-case name of AQ in string.
};

#endif //BART_SRC_AQDATA_AQ_BASE_H__
