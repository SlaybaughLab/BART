#include <deal.II/base/numbers.h>

#include <boost/algorithm/string.hpp>

#include <sstream>
#include <utility>
#include <iomanip>
#include <fstream>

#include "../include/problem_definition.h"

using namespace dealii;

template <int dim>
ProblemDefinition<dim>::ProblemDefinition (ParameterHandler &prm)
:
pi(numbers::PI),
transport_model_name(prm.get("transport model")),
discretization(prm.get("spatial discretization")),
n_azi(prm.get_integer("angular quadrature order")),
n_group(prm.get_integer("number of groups")),
n_material(prm.get_integer("number of materials")),
is_eigen_problem(prm.get_bool("do eigenvalue calculations")),
do_nda(prm.get_bool("do NDA")),
have_reflective_bc(prm.get_bool("have reflective BC")),
p_order(prm.get_integer("finite element polynomial degree")),
global_refinements(prm.get_integer("uniform refinements")),
output_namebase(prm.get("output file name base"))
{
  this->process_input (prm);
}

template <int dim>
ProblemDefinition<dim>::~ProblemDefinition()
{
}

template <int dim>
void ProblemDefinition<dim>::declare_parameters (ParameterHandler &prm)
{
  // our final strategy is to declare all possible entries
  // and then ignore some of them suggested by Wolfgang Bangerth
  // from Colorado State on 05-10-2017
  // The following are the basic parameters we need to define a problem
  {
    prm.declare_entry ("transport model", "ep", Patterns::Anything(), "valid names such as ep");
    prm.declare_entry ("angular quadrature order", "4", Patterns::Integer (), "Gauss-Chebyshev level-symmetric-like quadrature");
    prm.declare_entry ("number of groups", "1", Patterns::Integer (), "Number of groups in MG calculations");
    prm.declare_entry ("spatial discretization", "cfem", Patterns::Selection("DFEM|dfem|CFEM|cfem|dg|cg|DG|CG"), "USE DG or CG for spatial discretization");
    prm.declare_entry ("do eigenvalue calculations", "false", Patterns::Bool(), "Boolean to determine problem type");
    prm.declare_entry ("do NDA", "false", Patterns::Bool(), "Boolean to determine NDA or not");
    prm.declare_entry ("have reflective BC", "false", Patterns::Bool(), "");
    prm.declare_entry ("reflective boundary names", "", Patterns::List (Patterns::Anything ()), "must be lower cases of xmin,xmax,ymin,ymax,zmin,zmax");
    prm.declare_entry ("finite element polynomial degree", "1", Patterns::Integer(), "polynomial degree p for finite element");
    prm.declare_entry ("uniform refinements", "0", Patterns::Integer(), "number of uniform refinements desired");
    prm.declare_entry ("x, y, z max values of boundary locations", "", Patterns::List (Patterns::Double ()), "xmax, ymax, zmax of the boundaries, mins are zero");
    prm.declare_entry ("number of cells for x, y, z directions", "", Patterns::List (Patterns::Integer ()), "Geotry is hyper rectangle defined by how many cells exist per direction");
    prm.declare_entry ("number of materials", "1", Patterns::Integer (), "must be a positive integer");
    prm.declare_entry ("do print angular quadrature info", "true", Patterns::Bool(), "Boolean to determine if printing angular quadrature information");
    prm.declare_entry ("use explicit reflective boundary condition or not", "true", Patterns::Bool(), "");
    prm.declare_entry ("output file name base", "solu", Patterns::Anything(), "name base of the output file");
    // prm.declare_entry("material ID map", "", Patterns::List (Patterns::Integer ()), "Give material IDs for all blocks");
  }
  // FixIt: for current deal.II code, we don't consider reading mesh
  
  // Explanation: we brute-forcely declare as many entries as possible without read-in problem-definition
  // parameters. nmat and ngrp should both be large enough s.t. when reading starts, the real setting will
  // have entry-declaration
  prm.enter_subsection ("material ID map");
  {
    for (unsigned int z=0; z<z_levels; ++z)
      for (unsigned int y=0; y<y_levels; ++y)
      {
        std::ostringstream os;
        os << "cz " << z + 1 << ", cy " << y + 1;
        prm.declare_entry (os.str(), "", Patterns::List(Patterns::Integer()), "material IDs for a specific z level, y row");
      }
  }
  prm.leave_subsection ();
  
  prm.enter_subsection ("sigma_t, group=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm.declare_entry (os.str (), "", Patterns::List (Patterns::Double ()), "");
    }
  }
  prm.leave_subsection ();
  
  for (unsigned int m=0; m<nmat; ++m)
  {
    std::ostringstream os;
    os << "sigma_s, material " << m + 1;
    prm.enter_subsection (os.str());
    {
      for (unsigned int gin=0; gin<ngrp; ++gin)
      {
        std::ostringstream osm;
        osm << "g_in=" << gin + 1;
        prm.declare_entry (osm.str(), "", Patterns::List(Patterns::Double()), "multigroup sigma_s");
      }
    }
    prm.leave_subsection ();
  }
  
  prm.enter_subsection ("one-group sigma_t");
  {
    prm.declare_entry ("values", "", Patterns::List(Patterns::Double()), "");
  }
  prm.leave_subsection ();
  
  prm.enter_subsection ("one-group sigma_s");
  {
    prm.declare_entry ("values", "", Patterns::List(Patterns::Double()), "");
  }
  prm.leave_subsection ();
  
  prm.enter_subsection ("one-group Q");
  {
    prm.declare_entry ("values", "", Patterns::List(Patterns::Double()), "");
  }
  prm.leave_subsection ();
  
  prm.enter_subsection ("Q, group=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm.declare_entry (os.str (), "", Patterns::List (Patterns::Double ()), "");
    }
  }
  prm.leave_subsection ();
  
  // the following is for eigen problems
  prm.enter_subsection ("Fissile material IDs");
  {
    prm.declare_entry ("fissile material IDs", "", Patterns::List (Patterns::Integer ()), "");
  }
  prm.leave_subsection ();
  
  prm.enter_subsection ("ksi, group=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm.declare_entry(os.str(), "", Patterns::List(Patterns::Double()), "");
    }
  }
  prm.leave_subsection ();
  
  prm.enter_subsection ("nu_sigf, group=1 to G");
  {
    for (unsigned int m=0; m<nmat; ++m)
    {
      std::ostringstream os;
      os << "material " << m + 1;
      prm.declare_entry(os.str(), "", Patterns::List(Patterns::Double()), "");
    }
  }
  prm.leave_subsection ();
}

template <int dim>
std::string ProblemDefinition<dim>::get_transport_model ()
{
  return transport_model_name;
}

template <int dim>
std::string ProblemDefinition<dim>::get_transport_model (ParameterHandler &prm)
{
  std::string model_name = prm.get("transport model");
  return model_name;
}

template <int dim>
std::string ProblemDefinition<dim>::get_output_namebase ()
{
  return output_namebase;
}

template <int dim>
void ProblemDefinition<dim>::process_input (ParameterHandler &prm)
{
  preprocess_reflective_bc (prm);
  process_coordinate_information (prm);
  initialize_relative_position_to_id_map (prm);
  process_material_properties (prm);//???????????????????
  initialize_ref_bc_index ();
  produce_angular_quad ();
  initialize_component_index ();
}

template <int dim>
void ProblemDefinition<dim>::preprocess_reflective_bc (ParameterHandler &prm)
{
  if (have_reflective_bc)
  {
    is_explicit_reflective = prm.get_bool ("use explicit reflective boundary condition or not");
    // sorry, we need c++11
    std::map<std::string, unsigned int> bd_names_to_id {{"xmin",0},
      {"xmax",1},
      {"ymin",2},
      {"ymax",3},
      {"zmin",4},
      {"zmax",5}};
    std::vector<std::string> strings = Utilities::split_string_list (prm.get ("reflective boundary names"));
    AssertThrow (strings.size()>0,
                 ExcMessage("reflective boundary names have to be entered"));
    std::set<int> tmp;
    for (unsigned int i=0; i<strings.size (); ++i)
    {
      AssertThrow(bd_names_to_id.find(strings[i])!=bd_names_to_id.end(),
                  ExcMessage("Invalid reflective boundary name: use xmin, xmax, etc."));
      tmp.insert (bd_names_to_id[strings[i]]);
    }
    auto it = tmp.begin ();
    std::ostringstream os;
    os << "No valid reflective boundary name for " << dim << "D";
    AssertThrow(*it<2*dim,
                ExcMessage(os.str()));
    for (unsigned int i=0; i<2*dim; ++i)
    {
      if (tmp.count (i))
        is_reflective_bc[i] = true;
      else
        is_reflective_bc[i] = false;
    }
  }
}

template <int dim>
bool ProblemDefinition<dim>::get_explicit_reflective_bool ()
{
  return is_explicit_reflective;
}

template <int dim>
std::unordered_map<unsigned int, bool>
ProblemDefinition<dim>::get_reflective_bc_map ()
{
  return is_reflective_bc;
}

template <int dim>
bool ProblemDefinition<dim>::get_print_sn_quad_bool ()
{
  return do_print_sn_quad;
}

template <int dim>
void ProblemDefinition<dim>::process_coordinate_information (ParameterHandler &prm)
{
  // max values for all axis
  {
    std::vector<std::string> strings = Utilities::split_string_list (prm.get ("x, y, z max values of boundary locations"));
    AssertThrow (strings.size()>=dim,
                 ExcMessage("Number of axis max values must be the same as dimension"));
    for (unsigned int i=0; i<dim; ++i)
      axis_max_values.push_back (std::atof (strings[i].c_str()));
  }
  
  // read in number of cells and get cell sizes along axes
  {
    std::vector<std::string> strings = Utilities::split_string_list (prm.get ("number of cells for x, y, z directions"));
    AssertThrow (strings.size()>=dim,
                 ExcMessage ("Entries for numbers of cells should be equal to dimension"));
    std::vector<unsigned int> cells_per_dir;
    std::vector<std::vector<double> > spacings;
    for (unsigned int d=0; d<dim; ++d)
    {
      ncell_per_dir.push_back (std::atoi (strings[d].c_str ()));
      cell_size_all_dir.push_back (axis_max_values[d]/ncell_per_dir[d]);
    }
  }
}

template <int dim>
std::vector<double>
ProblemDefinition<dim>::get_axis_maxes ()
{
  return axis_max_values;
}

template <int dim>
std::vector<unsigned int>
ProblemDefinition<dim>::get_ncells ()
{
  return ncell_per_dir;
}

template <int dim>
std::vector<double>
ProblemDefinition<dim>::get_cell_sizes ()
{
  return cell_size_all_dir;
}

template <int dim>
std::map<std::vector<unsigned int>, unsigned int>
ProblemDefinition<dim>::get_id_map ()
{
  return relative_position_to_id;
}

template <int dim>
void ProblemDefinition<dim>::initialize_relative_position_to_id_map (ParameterHandler &prm)
{
  prm.enter_subsection ("material ID map");
  {
    unsigned int ncell_z = dim==3?ncell_per_dir[2]:1;
    unsigned int ncell_y = dim>=2?ncell_per_dir[1]:1;
    for (unsigned int z=0; z<ncell_z; ++z)
      for (unsigned int y=0; y<ncell_y; ++y)
      {
        std::ostringstream os;
        os << "cz " << z + 1 << ", cy " << y + 1;
        std::vector<std::string> strings = Utilities::split_string_list (prm.get(os.str()));
        for (unsigned int x=0; x<ncell_per_dir[0]; ++x)
        {
          std::vector<unsigned int> tmp {x, y, z};
          std::cout << x << ", " << y << ", " << z << ", "
          << std::atoi (strings[x].c_str()) - 1 << ", is the id" << std::endl;
          AssertThrow(std::atoi(strings[x].c_str())>0,
                      ExcMessage("material ID must be larger than 0"));
          relative_position_to_id[tmp] = std::atoi (strings[x].c_str()) - 1;
        }
      }
    std::vector<unsigned int> t1{0,0,0}, t2{0,1,0}, t3{1,0,0}, t4{1,1,0};
    std::cout << "tests: " << 0 << ", " << 0 << ", " << relative_position_to_id[t1] << std::endl;
    std::cout << "tests: " << 0 << ", " << 1 << ", " << relative_position_to_id[t2] << std::endl;
    std::cout << "tests: " << 1 << ", " << 0 << ", " << relative_position_to_id[t3] << std::endl;
    std::cout << "tests: " << 1 << ", " << 1 << ", " << relative_position_to_id[t4] << std::endl;
  }
  prm.leave_subsection ();
}

template <int dim>
void ProblemDefinition<dim>::process_material_properties (ParameterHandler &prm)
{
  initialize_relative_position_to_id_map (prm);
  if (n_group>1)
  {
    // This block takes in sigts
    prm.enter_subsection ("sigma_t, group=1 to G");
    {
      for (unsigned int m=0; m<n_material; ++m)
      {
        std::ostringstream os;
        os << "material " << m + 1;
        std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
        AssertThrow (strings.size () == n_group,
                     ExcMessage ("n_group is not equal to group number of sigma_t"));
        std::vector<double> tmp, inv_tmp;
        for (unsigned int g=0; g<n_group; ++g)
        {
          tmp.push_back (std::atof (strings[g].c_str ()));
          inv_tmp.push_back (1.0/tmp[g]);
        }
        all_sigt.push_back (tmp);
        all_inv_sigt.push_back (inv_tmp);
      }
    }
    prm.leave_subsection ();
    
    // This block takes in scattering transfer matrices
    all_sigs.resize (n_material);
    all_sigs_per_ster.resize (n_material);
    for (unsigned int m=0; m<n_material; ++m)
    {
      std::ostringstream osm;
      osm << "sigma_s, material " << m + 1;
      std::vector<std::vector<double> >  tmp_sigs (n_group, std::vector<double>(n_group));
      std::vector<std::vector<double> >  tmp_sigs_per_ster (n_group, std::vector<double>(n_group));
      prm.enter_subsection (osm.str());
      {
        for (unsigned int gin=0; gin<n_group; ++gin)
        {
          std::ostringstream os;
          os << "g_in=" << gin + 1;
          std::vector<std::string> strings = Utilities::split_string_list (prm.get(os.str()));
          AssertThrow (strings.size()==n_group,
                       ExcMessage("sigma_s should have n_group entries per in group"));
          for (unsigned int g=0; g<n_group; ++g)
          {
            tmp_sigs[gin][g] = std::atof (strings[g].c_str());
            tmp_sigs_per_ster[gin][g] = std::atof (strings[g].c_str()) / (4.0 * pi);
          }
        }
      }
      prm.leave_subsection ();
      all_sigs[m] = tmp_sigs;
      all_sigs_per_ster[m] = tmp_sigs_per_ster;
    }
    
    if (!is_eigen_problem)
    {
      prm.enter_subsection ("Q, group=1 to G");
      {
        for (unsigned int m=0; m<n_material; ++m)
        {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
          AssertThrow (strings.size () == n_group,
                       ExcMessage ("n_group is not equal to group number of Q"));
          std::vector<double> tmp_q, tmp_q_per_ster;
          for (unsigned int g=0; g<n_group; ++g)
          {
            tmp_q.push_back (std::atof (strings[g].c_str ()));
            tmp_q_per_ster.push_back (tmp_q[g] / (4.0 * pi));
          }
          all_q.push_back (tmp_q);
          all_q_per_ster.push_back (tmp_q_per_ster);
        }
      }
      prm.leave_subsection ();
    }
  }
  else
  {
    prm.enter_subsection ("one-group sigma_t");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()>=n_material,
                   ExcMessage("One-group sigma_t should have N_material entries"));
      
      for (unsigned int m=0; m<n_material; ++m)
      {
        // sorry, c++11 again.
        std::vector<double> tmp = {std::atof (strings[m].c_str())};
        std::vector<double> inv_tmp = {1.0 / std::atof (strings[m].c_str())};
        
        all_sigt.push_back (tmp);
        all_inv_sigt.push_back (inv_tmp);
      }
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("one-group sigma_s");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group sigma_s should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        // sorry, c++11 again.
        std::vector<std::vector<double> > tmp {std::vector<double> {std::atof(strings[m].c_str())}};
        all_sigs.push_back (tmp);
        std::vector<std::vector<double> > tmp_per_ster {std::vector<double> {std::atof(strings[m].c_str()) / (4.0*pi)}};
        all_sigs_per_ster.push_back (tmp_per_ster);
      }
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("one-group Q");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group Q should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        std::vector<double> tmp = {std::atof (strings[m].c_str())};
        all_q.push_back (tmp);
        std::vector<double> tmp_per_ster = {std::atof (strings[m].c_str()) / (4.0 * pi)};
        all_q_per_ster.push_back (tmp_per_ster);
      }
    }
    prm.leave_subsection ();
  }
  
  // This block is for eigenvalue problems
  if (is_eigen_problem)
  {
    process_eigen_material_properties (prm);
  }
}

template <int dim>
void ProblemDefinition<dim>::process_eigen_material_properties (ParameterHandler &prm)
{
  prm.enter_subsection ("Fissile material IDs");
  {
    std::ostringstream os;
    os << "fissile material IDs";
    std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
    AssertThrow (strings.size () > 0,
                 ExcMessage ("Fissile material IDs must be inserted for eigen problems"));
    // std::set<int> fissile_ids;
    for (unsigned int i=0; i<strings.size(); ++i)
      fissile_ids.insert (std::atoi (strings[i].c_str ()));
  }
  prm.leave_subsection ();
  
  for (unsigned int m=0; m<n_material; ++m)
  {
    if (fissile_ids.count(m))
      is_material_fissile[m] = true;
    else
      is_material_fissile[m] = false;
  }
  AssertThrow (!is_material_fissile.empty (),
               ExcMessage ("Please specify at least one valid ID for fissile materials"));
  
  if (n_group>1)
  {
    prm.enter_subsection ("ksi, group=1 to G");
    {
      for (unsigned int m=0; m<n_material;++m)
      {
        std::vector<double> tmp(n_group,0.0);
        if (is_material_fissile[m])
        {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
          AssertThrow (strings.size () == n_group,
                       ExcMessage ("n_group is not equal to group number of ksi"));
          std::vector<double> tmp;
          for (unsigned int g=0; g<n_group; ++g)
            tmp[g] = std::atof (strings[g].c_str ());
        }
        all_ksi.push_back (tmp);
      }
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("nu_sigf, group=1 to G");
    {
      for (unsigned int m=0; m<n_material;++m)
      {
        std::vector<double> tmp(n_group,0.0);
        if (is_material_fissile[m])
        {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
          AssertThrow (strings.size () == n_group,
                       ExcMessage ("n_group is not equal to group number of nusigf"));
          std::vector<double> tmp;
          for (unsigned int g=0; g<n_group; ++g)
            tmp[g] = std::atof (strings[g].c_str ());
        }
        all_nusigf.push_back (tmp);
      }
    }
    prm.leave_subsection ();
  }
  else
  {
    prm.enter_subsection ("one-group ksi");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group ksi should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        // sorry, c++11 again.
        std::vector<double> tmp {is_material_fissile[m] ? std::atof (strings[m].c_str()) : 0.0};
        all_sigt.push_back (tmp);
      }
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("one-group nu_sigf");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group sigma_t should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        std::vector<double> tmp {is_material_fissile[m] ? std::atof (strings[m].c_str()) : 0.0};
        all_sigt.push_back (tmp);
      }
    }
    prm.leave_subsection ();
  }
  
  for (unsigned int m=0; m<n_material; ++m)
  {
    std::vector<std::vector<double> >  tmp (n_group, std::vector<double>(n_group));
    std::vector<std::vector<double> >  tmp_per_ster (n_group, std::vector<double>(n_group));
    if (is_material_fissile[m])
      for (unsigned int gin=0; gin<n_group; ++gin)
        for (unsigned int g=0; g<n_group; ++g)
        {
          tmp[gin][g] = all_ksi[m][g] * all_nusigf[m][gin];
          tmp_per_ster[gin][g] = tmp[gin][g] / (4.0 * pi);
        }
    all_ksi_nusigf.push_back (tmp_per_ster);
    all_ksi_nusigf_per_ster.push_back (tmp_per_ster);
  }
}

template <int dim>
std::unordered_map<unsigned int, bool> ProblemDefinition<dim>::get_fissile_id_map ()
{
  return is_material_fissile;
}

template <int dim>
void ProblemDefinition<dim>::initialize_ref_bc_index ()
{
  if (have_reflective_bc)
  {
    AssertThrow (dim>1,
                 ExcMessage("1D cases are not implemented for now."));
    std::vector<Tensor<1, dim> > boundary_normal_vectors;
    boundary_normal_vectors.resize (2*dim);
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
        Tensor<1, dim> out_angle = (omega_i[i_dir]
                                    *
                                    (1.0 - 2.0 * (boundary_normal_vectors[i] * omega_i[i_dir])));
        for (unsigned int r_dir=0; r_dir<n_dir; ++r_dir)
        {
          Tensor<1, dim> d_dir = out_angle;
          Tensor<1, dim> d_minus_dir = out_angle;
          d_dir -= omega_i[r_dir];
          d_minus_dir += omega_i[r_dir];
          if (d_dir.norm ()<1.0e-13 || d_minus_dir.norm ()<1.0e-13)
            reflective_direction_index.insert (std::make_pair (std::make_pair (i, i_dir), r_dir));
        }
      }
  }
}

template <int dim>
std::map<std::pair<unsigned int, unsigned int>, unsigned int>
ProblemDefinition<dim>::get_reflective_direction_index_map ()
{
  return reflective_direction_index;
}

template <int dim>
void ProblemDefinition<dim>::produce_angular_quad ()
{
  Assert (n_azi%2==0, ExcDimensionMismatch(n_azi%2, 1));
  
  QGauss<1> mu_quad (n_azi);
  
  
  if (dim==1)
  {
    total_angle = 2.0;
    for (unsigned int i=0; i<n_azi/2; ++i)
    {
      Tensor <1, dim> tmp;
      tmp[0] = mu_quad.point(i)[0] * 2.0 - 1.0;
      omega_i.push_back (tmp);
      wi.push_back (mu_quad.weight(i));
    }
    n_dir = n_azi / 2;
  }
  
  if (dim==2)
  {
    total_angle = 4.0 * pi;
    double aha=0;
    for (unsigned int i=n_azi/2; i<n_azi; ++i)
    {
      Tensor<1, dim> tmp;
      double mut = mu_quad.point(i)[0]*2.0 - 1.0;
      unsigned int level_angle_num = 4 * (n_azi - i);
      double delta_level_angle = 2.0 * pi / level_angle_num;
      double level_weight = mu_quad.weight(i) * 16.0 * pi;
      
      for (unsigned int j=0; j<level_angle_num; ++j)
      {
        double angle = 0.5 * delta_level_angle + (double)(j) * delta_level_angle;
        tmp[0] = std::sqrt (1.0 - mut * mut) * cos (angle);
        tmp[1] = std::sqrt (1.0 - mut * mut) * sin (angle);
        // solely for EP quadratures in 2D.
        if (tmp[0]>0)
        {
          omega_i.push_back(tmp);
          double point_wt = level_weight / level_angle_num;
          wi.push_back(point_wt);
          Tensor<1,3> tmp2;
          tmp2[0] = tmp[0];
          tmp2[1] = tmp[1];
          tmp2[2] = mut;
          omega_with_mu.push_back (tmp2);
        }
        
      }
      
    }
    n_dir = n_azi * (n_azi + 2) / 4;
  }
  
  if (dim==3)
  {
    total_angle = 4.0 * pi;
    
    for (unsigned int i=0; i<n_azi/2; ++i)
    {
      unsigned int level_angle_num = i < n_azi / 2 ? 4 * (i + 1) : 4 * (n_azi - i);
      double delta_level_angle = 2.0 * pi / level_angle_num;
      double level_weight = mu_quad.weight(i) * 8.0 * pi;
      Tensor<1, dim> tmp;
      double mut = mu_quad.point(i)[0] * 2.0 - 1.0;
      
      for (unsigned int j=0; j<level_angle_num; ++j)
      {
        double angle = 0.5 * delta_level_angle + (double)(j) * delta_level_angle;
        tmp[0] = std::sqrt(1.0 - mut * mut) * cos(angle);
        tmp[1] = std::sqrt(1.0 - mut * mut) * sin(angle);
        tmp[2] = mut;
        omega_i.push_back(tmp);
        double point_wt = level_weight / level_angle_num;
        wi.push_back(point_wt);
        Tensor<1,3> tmp2;
        tmp2[0] = tmp[0];
        tmp2[1] = tmp[1];
        tmp2[2] = mut;
        omega_with_mu.push_back (tmp2);
      }
      
    }
    n_dir = n_azi * (n_azi + 2) / 2;
  }
  AssertThrow (n_dir==wi.size(),
               ExcMessage("calculated number of angles should be the same as number of angular weights"));
  AssertThrow (n_dir==omega_i.size(),
               ExcMessage("calculated number of angles should be the same as number of angles"));
  
  n_total_ho_vars = n_dir * n_group;
  
  /*
   for (unsigned int i=0; i<n_total_ho_vars; ++i)
   {
   FEValuesExtractors::Scalar tmp(i);
   comp.push_back(tmp);
   }
   */
  
  // estimate tensor norm to do penalty method
  for (unsigned int i=0; i<n_dir; ++i)
  {
    Tensor<2, dim> tensor_tmp = outer_product(omega_i[i], omega_i[i]);
    //double norm = 10.;
    //double norm = l1_norm (tensor_tmp);
    double norm = tensor_tmp.norm();
    //double norm = linfty_norm (tensor_tmp);
    std::cout << "tensor norm " << norm << std::endl;
    tensor_norms.push_back(norm);
  }
}

template <int dim>
unsigned int ProblemDefinition<dim>::get_sn_order ()
{
  return n_azi;
}

template <int dim>
unsigned int ProblemDefinition<dim>::get_n_dir ()
{
  return n_dir;
}

template <int dim>
unsigned int ProblemDefinition<dim>::get_n_group ()
{
  return n_group;
}

template <int dim>
unsigned int ProblemDefinition<dim>::get_n_total_ho_vars ()
{
  return n_total_ho_vars;
}

template <int dim>
std::vector<double> ProblemDefinition<dim>::get_angular_weights ()
{
  return wi;
}

template <int dim>
std::vector<Tensor<1, dim> > ProblemDefinition<dim>::get_all_directions ()
{
  return omega_i;
}

template <int dim>
std::vector<double> ProblemDefinition<dim>::get_tensor_norms ()
{
  return tensor_norms;
}

template <int dim>
bool ProblemDefinition<dim>::get_nda_bool ()
{
  return do_nda;
}

template <int dim>
bool ProblemDefinition<dim>::get_eigen_problem_bool ()
{
  return is_eigen_problem;
}

template <int dim>
std::vector<std::vector<double> > ProblemDefinition<dim>::get_sigma_t ()
{
  return all_sigt;
}

template <int dim>
std::vector<std::vector<double> > ProblemDefinition<dim>::get_inv_sigma_t ()
{
  return all_inv_sigt;
}

template <int dim>
std::vector<std::vector<double> > ProblemDefinition<dim>::get_q ()
{
  return all_q;
}

template <int dim>
std::vector<std::vector<double> > ProblemDefinition<dim>::get_q_per_ster ()
{
  return all_q_per_ster;
}

template <int dim>
std::vector<std::vector<std::vector<double> > > ProblemDefinition<dim>::get_sigma_s ()
{
  return all_sigs;
}

template <int dim>
std::vector<std::vector<std::vector<double> > > ProblemDefinition<dim>::get_sigma_s_per_ster ()
{
  return all_sigs_per_ster;
}

template <int dim>
std::vector<std::vector<std::vector<double> > > ProblemDefinition<dim>::get_ksi_nusigf ()
{
  return all_ksi_nusigf;
}

template <int dim>
std::vector<std::vector<std::vector<double> > > ProblemDefinition<dim>::get_ksi_nusigf_per_ster ()
{
  return all_ksi_nusigf_per_ster;
}

template <int dim>
std::vector<std::vector<double> > ProblemDefinition<dim>::get_nusigf ()
{
  return all_nusigf;
}

template <int dim>
bool ProblemDefinition<dim>::get_reflective_bool ()
{
  return have_reflective_bc;
}

template <int dim>
unsigned int ProblemDefinition<dim>::get_fe_order ()
{
  return p_order;
}

template <int dim>
unsigned int ProblemDefinition<dim>::get_uniform_refinement ()
{
  return global_refinements;
}

template <int dim>
void ProblemDefinition<dim>::initialize_component_index ()
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
std::map<std::pair<unsigned int, unsigned int>, unsigned int>
ProblemDefinition<dim>::get_component_index_map ()
{
  return component_index;
}

template <int dim>
std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> >
ProblemDefinition<dim>::get_inv_component_map ()
{
  return inverse_component_index;
}

template <int dim>
std::string ProblemDefinition<dim>::get_discretization ()
{
  return discretization;
}

template <int dim>
void ProblemDefinition<dim>::print_angular_quad ()
{
  if (do_print_sn_quad)
  {
    if (dim==1)
      AssertThrow (dim>=2,
                   ExcMessage("1D is not implemented yet."));
    
    std::ofstream quadr;
    quadr.open("quadr.txt");
    quadr << "Dim = " << dim << ", SN order = " << n_azi << std::endl;
    quadr << std::setfill ('+') << std::setw (57) << std::endl;
    quadr << std::setfill ('+') << std::setw (57) << std::endl;
    quadr << "Weights | Omega_z | Omega_x | Omega_y" << std::endl;
    quadr << std::setfill ('-') << std::setw (57) << std::endl;
    quadr << std::setfill ('-') << std::setw (57) << std::endl;
    for (unsigned int i=0; i<omega_i.size(); ++i)
    {
      quadr << std::fixed << std::setprecision (15);
      quadr << wi[i] << ", ";
      std::cout << wi[i] << std::endl;
      quadr << omega_with_mu[i][2] << ", ";
      quadr << omega_with_mu[i][0] << ", ";
      quadr << omega_with_mu[i][1] << std::endl;
      quadr << std::setfill ('-') << std::setw (57) << std::endl;
    }
    quadr.close ();
  }
}

// explicit instantiation to avoid linking error
template class ProblemDefinition<2>;
template class ProblemDefinition<3>;
