#include <deal.II/base/numbers.h>

#include <fstream>

#include "../include/aq_base.h"

template <int dim>
AQBase<dim>::AQBase (ParameterHandler &prm)
:
pi(numbers::PI),
transport_model_name(prm.get("transport model")),
aq_name(prm.get("angular quadrature name")),
n_group(prm.get_integer("number of groups")),
n_azi(prm.get_integer("angular quadrature order")),
have_reflective_bc(prm.get_bool("have reflective BC"))
{
}

template <int dim>
AQBase<dim>::~AQBase ()
{
}

template <int dim>
void AQBase<dim>::make_aq (ParameterHandler &prm)
{
  if (transport_model_name=="ep")
    discretization = prm.get ("spatial discretization");
  produce_angular_quad ();
  initialize_component_index ();
  initialize_ref_bc_index ();
}

template <int dim>
void AQBase<dim>::initialize_ref_bc_index ()
{
  // Note: here we assume square domain and assume user either
  // uses deal.II generated mesh or mesh from gmsh with proper
  // boundary IDs setup: {0,1,2,3,4,5} for {xmin,xmax,ymin,ymax,zmin,zmax}
  if (have_reflective_bc)
  {
    AssertThrow (dim>1,
                 ExcMessage("1D cases are not implemented for now."));
    std::vector<Tensor<1, dim> > boundary_normal_vectors;
    boundary_normal_vectors.resize (2*dim);
    // All boundary normal vectors are assume to be parallel to axes
    // Then, only one component in each normal vector is nonzero
    if (dim==2)
    {
      boundary_normal_vectors[0][0] = -1.0;
      boundary_normal_vectors[1][0] = 1.0;
      boundary_normal_vectors[2][1] = -1.0;
      boundary_normal_vectors[3][1] = 1.0;
    }
    if (dim==3)
    {
      boundary_normal_vectors[0][0] = -1.0;
      boundary_normal_vectors[1][0] = 1.0;
      boundary_normal_vectors[2][1] = -1.0;
      boundary_normal_vectors[3][1] = 1.0;
      boundary_normal_vectors[4][2] = -1.0;
      boundary_normal_vectors[5][2] = 1.0;
    }
    for (unsigned int i=0; i<2*dim; ++i)
      for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
      {
        Tensor<1, dim> out_angle =
        (omega_i[i_dir] *
         (1.0 - 2.0 * (boundary_normal_vectors[i] * omega_i[i_dir])));
        for (unsigned int r_dir=0; r_dir<n_dir; ++r_dir)
        {
          Tensor<1, dim> d_dir = out_angle;
          Tensor<1, dim> d_minus_dir = out_angle;
          d_dir -= omega_i[r_dir];
          d_minus_dir += omega_i[r_dir];
          if (transport_model_name=="ep")
            if (d_dir.norm ()<1.0e-13 || d_minus_dir.norm ()<1.0e-13)
              reflective_direction_index.insert (std::make_pair (std::make_pair (i, i_dir), r_dir));
            else
              if (d_dir.norm ()<1.0e-13)
                reflective_direction_index.insert (std::make_pair (std::make_pair (i, i_dir), r_dir));
        }
      }
  }
}

template <int dim>
void AQBase<dim>::produce_angular_quad ()
{
}

template <int dim>
void AQBase<dim>::initialize_component_index ()
{
  // initialize the map from (direction, group) to component indices
  unsigned int ind = 0;
  for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    for (unsigned int g=0; g<n_group; ++g)
    {
      std::pair<unsigned int, unsigned int> key (i_dir, g);
      component_index[key] = ind;
      inverse_component_index[ind] = key;
      ind += 1;
    }
}

template <int dim>
void AQBase<dim>::print_angular_quad ()
{
  AssertThrow (dim>=2,
               ExcMessage("1D is not implemented yet."));
  std::ofstream quadr;
  quadr.open("aq.txt");
  quadr << "transport model: " << transport_model_name
  << "; quadrature name: " << produce_quadrature_name () << std::endl;
  quadr << "Dim = " << dim << ", SN order = " << n_azi << std::endl;
  quadr << "Weights | Omega_x | Omega_y | mu" << std::endl;
  for (unsigned int i=0; i<omega_i.size(); ++i)
  {
    double mu = std::sqrt (1-(omega_i[i][0]*omega_i[i][0]+omega_i[i][1]*omega_i[i][1]));
    quadr << std::fixed << std::setprecision (15);
    quadr << wi[i] << ", ";
    quadr << omega_i[i][0] << ", ";
    quadr << omega_i[i][1] << ", ";
    quadr << mu << std::endl;
  }
  quadr.close ();
}

template <int dim>
std::string AQBase<dim>::produce_quadrature_name ()
{
  AssertThrow (aq_name.size()>0,
               ExcMessage("aq name has to be assigned before getting a physical name"));
  // ToDo: more quadrature name producers
  if (aq_name=="lsgc")
    return "Level Symmetric Gauss Chebyshev";
}

//public member functions to retrieve private and protected variables
template <int dim>
std::map<std::pair<unsigned int, unsigned int>, unsigned int>
AQBase<dim>::get_component_index_map ()
{
  return component_index;
}

template <int dim>
std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> >
AQBase<dim>::get_inv_component_map ()
{
  return inverse_component_index;
}

template <int dim>
std::map<std::pair<unsigned int, unsigned int>, unsigned int>
AQBase<dim>::get_reflective_direction_index_map ()
{
  return reflective_direction_index;
}

template <int dim>
unsigned int AQBase<dim>::get_sn_order ()
{
  return n_azi;
}

template <int dim>
unsigned int AQBase<dim>::get_n_dir ()
{
  return n_dir;
}

template <int dim>
unsigned int AQBase<dim>::get_n_total_ho_vars ()
{
  return n_total_ho_vars;
}

template <int dim>
std::vector<double> AQBase<dim>::get_angular_weights ()
{
  return wi;
}

template <int dim>
std::vector<Tensor<1, dim> > AQBase<dim>::get_all_directions ()
{
  return omega_i;
}

template <int dim>
std::vector<double> AQBase<dim>::get_tensor_norms ()
{
  return tensor_norms;
}

template class AQBase<2>;
template class AQBase<3>;
