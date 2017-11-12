#include <deal.II/base/numbers.h>
#include <deal.II/base/utilities.h>

#include "material_properties.h"

MaterialProperties::MaterialProperties (ParameterHandler &prm)
:
pi(numbers::PI),
is_eigen_problem(prm.get_bool("do eigenvalue calculations")),
do_nda(prm.get_bool("do NDA")),
n_group(prm.get_integer("number of groups")),
n_material(prm.get_integer("number of materials"))
{
  process_material_properties (prm);
}

MaterialProperties::~MaterialProperties ()
{
}

void MaterialProperties::process_material_properties
(ParameterHandler &prm)
{
  if (n_group>1)
  {
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
  }// MG
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
        std::vector<std::vector<double> > tmp {std::vector<double> {std::atof(strings[m].c_str())}};
        all_sigs.push_back (tmp);
        std::vector<std::vector<double> > tmp_per_ster {std::vector<double> {std::atof(strings[m].c_str()) / (4.0*pi)}};
        all_sigs_per_ster.push_back (tmp_per_ster);
      }
    }
    prm.leave_subsection ();
    
    if (!is_eigen_problem)
    {
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
  }
  
  // This block is for eigenvalue problems
  if (is_eigen_problem)
  {
    process_eigen_material_properties (prm);
  }
}

void MaterialProperties::process_eigen_material_properties
(ParameterHandler &prm)
{
  prm.enter_subsection ("fissile material IDs");
  {
    std::ostringstream os;
    os << "fissile material ids";
    std::vector<std::string> strings = Utilities::split_string_list (prm.get (os.str ()));
    AssertThrow (strings.size () > 0,
                 ExcMessage ("Fissile material IDs must be inserted for eigen problems"));
    // std::set<int> fissile_ids;
    for (unsigned int i=0; i<strings.size(); ++i)
      fissile_ids.insert (std::atoi(strings[i].c_str ())-1);
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
    prm.enter_subsection ("chi, group=1 to G");
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
                       ExcMessage ("n_group is not equal to group number of chi"));
          std::vector<double> tmp;
          for (unsigned int g=0; g<n_group; ++g)
            tmp[g] = std::atof (strings[g].c_str ());
        }
        all_chi.push_back (tmp);
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
    prm.enter_subsection ("one-group chi");
    {
      std::vector<std::string> strings = Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material,
                   ExcMessage("One-group chi should have N_material entries"));
      std::vector<double> tmp_sigt (n_material);
      for (unsigned int m=0; m<n_material; ++m)
      {
        // sorry, c++11 again.
        std::vector<double> tmp {is_material_fissile[m] ? std::atof (strings[m].c_str()) : 0.0};
        all_chi.push_back (tmp);
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
        std::vector<double> tmp
        {is_material_fissile[m] ? std::atof (strings[m].c_str()) : 0.0};
        all_nusigf.push_back (tmp);
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
          tmp[gin][g] = all_chi[m][g] * all_nusigf[m][gin];
          tmp_per_ster[gin][g] = tmp[gin][g] / (4.0 * pi);
        }
    all_chi_nusigf.push_back (tmp_per_ster);
    all_chi_nusigf_per_ster.push_back (tmp_per_ster);
  }
}

std::vector<std::vector<double> > MaterialProperties::get_sigma_t ()
{
  return all_sigt;
}

std::vector<std::vector<double> > MaterialProperties::get_inv_sigma_t ()
{
  return all_inv_sigt;
}

std::vector<std::vector<double> > MaterialProperties::get_q ()
{
  return all_q;
}

std::vector<std::vector<double> > MaterialProperties::get_q_per_ster ()
{
  return all_q_per_ster;
}

std::vector<std::vector<std::vector<double> > > MaterialProperties::get_sigma_s ()
{
  return all_sigs;
}

std::vector<std::vector<std::vector<double> > > MaterialProperties::get_sigma_s_per_ster ()
{
  return all_sigs_per_ster;
}

std::vector<std::vector<std::vector<double> > > MaterialProperties::get_chi_nusigf ()
{
  return all_chi_nusigf;
}

std::vector<std::vector<std::vector<double> > > MaterialProperties::get_chi_nusigf_per_ster ()
{
  return all_chi_nusigf_per_ster;
}

std::vector<std::vector<double> > MaterialProperties::get_nusigf ()
{
  return all_nusigf;
}

std::unordered_map<unsigned int, bool>
MaterialProperties::get_fissile_id_map ()
{
  return is_material_fissile;
}
