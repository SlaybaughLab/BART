#include "../common/numbers.h"
#include "materials.h"

#include <deal.II/base/utilities.h>

Materials::Materials (dealii::ParameterHandler &prm)
    : is_eigen_problem_(prm.get_bool("do eigenvalue calculations")),
      do_nda_(prm.get_bool("do nda")),
      n_group_(prm.get_integer("number of groups")),
      n_material_(prm.get_integer("number of materials")) {
  ProcessMaterials(prm);
}

Materials::~Materials () {}

void Materials::ProcessMaterials (dealii::ParameterHandler &prm) {
  if (n_group_>1) {
    prm.enter_subsection ("sigma_t, group=1 to G");
    {
      for (int m=0; m<n_material_; ++m) {
        std::ostringstream os;
        os << "material " << m + 1;
        std::vector<std::string> strings =
            dealii::Utilities::split_string_list (prm.get(os.str()));
        AssertThrow (strings.size () == n_group_,
            dealii::ExcMessage ("n_group_ is not equal to group number of sigma_t"));
        std::vector<double> tmp, inv_tmp;
        for (int g=0; g<n_group_; ++g) {
          tmp.push_back (std::atof (strings[g].c_str ()));
          inv_tmp.push_back (1.0/tmp[g]);
        }
        sigt_[m] = tmp;
        inv_sigt_[m] = inv_tmp;
      }
    }
    prm.leave_subsection ();

    // This block takes in scattering transfer matrices
    for (int m=0; m<n_material_; ++m) {
      std::ostringstream osm;
      osm << "sigma_s, material " << m + 1;
      dealii::FullMatrix<double> tmp_sigs (n_group_, n_group_);
      prm.enter_subsection (osm.str());
      {
        for (int gin=0; gin<n_group_; ++gin) {
          std::ostringstream os;
          os << "g_in=" << gin + 1;
          std::vector<std::string> strings =
              dealii::Utilities::split_string_list (prm.get(os.str()));
          AssertThrow (strings.size()==n_group_,
                       dealii::ExcMessage("sigma_s should have n_group_ entries per in group"));
          for (int g=0; g<n_group_; ++g) {
            tmp_sigs(gin, g) = std::atof (strings[g].c_str());
          }
        }
      }
      prm.leave_subsection ();
      sigs_[m] = tmp_sigs;
      sigs_per_ster_[m] = tmp_sigs;
      sigs_per_ster_[m] *= bconst::kInvFourPi;
    }

    if (!is_eigen_problem_) {
      prm.enter_subsection ("Q, group=1 to G");
      {
        for (int m=0; m<n_material_; ++m) {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings = dealii::Utilities::split_string_list (prm.get (os.str ()));
          AssertThrow (strings.size () == n_group_,
                       dealii::ExcMessage ("n_group_ is not equal to group number of Q"));
          std::vector<double> tmp_q;
          for (int g=0; g<n_group_; ++g)
            tmp_q.push_back (std::atof (strings[g].c_str ()));
          q_[m] = tmp_q;
          q_per_ster_[m] = tmp_q;
          std::for_each(q_per_ster_[m].begin(), q_per_ster_[m].end(),
              [&](double &val){return val*bconst::kInvFourPi;});
        }
      }
      prm.leave_subsection ();
    }
  } else {
    prm.enter_subsection ("one-group sigma_t");
    {
      std::vector<std::string> strings =
          dealii::Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()>=n_material_,
          dealii::ExcMessage("One-group sigma_t should have n_material_ entries"));

      for (int m=0; m<n_material_; ++m) {
        std::vector<double> tmp = {std::atof (strings[m].c_str())};
        std::vector<double> inv_tmp = {1.0 / std::atof (strings[m].c_str())};
        sigt_[m] = tmp;
        inv_sigt_[m] = inv_tmp;
      }
    }
    prm.leave_subsection ();

    prm.enter_subsection ("one-group sigma_s");
    {
      std::vector<std::string> strings =
          dealii::Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material_,
          dealii::ExcMessage("One-group sigma_s should have n_material_ entries"));
      for (int m=0; m<n_material_; ++m) {
        dealii::FullMatrix<double> tmp(1,1);
        tmp(0,0) = std::atof(strings[m].c_str());
        sigs_[m] = tmp;
        sigs_per_ster_[m] = tmp;
        sigs_per_ster_[m] *= bconst::kInvFourPi;
      }
    }
    prm.leave_subsection ();

    if (!is_eigen_problem_) {
      prm.enter_subsection ("one-group Q");
      {
        std::vector<std::string> strings = dealii::Utilities::split_string_list (prm.get("values"));
        AssertThrow (strings.size()==n_material_,
                     dealii::ExcMessage("One-group Q should have n_material_ entries"));
        std::vector<double> tmp_sigt (n_material_);
        for (int m=0; m<n_material_; ++m) {
          std::vector<double> tmp = {std::atof (strings[m].c_str())};
          q_[m] = tmp;
          std::vector<double> tmp_per_ster = {std::atof(strings[m].c_str())/bconst::kFourPi};
          q_per_ster_[m] = tmp_per_ster;
        }
      }
      prm.leave_subsection ();
    }
  }

  // This block is for eigenvalue problems
  if (is_eigen_problem_)
    ProcessEigenMaterials (prm);
}

void Materials::ProcessEigenMaterials (dealii::ParameterHandler &prm) {
  prm.enter_subsection ("fissile material IDs");
  {
    std::ostringstream os;
    os << "fissile material ids";
    std::vector<std::string> strings = dealii::Utilities::split_string_list (prm.get (os.str ()));
    AssertThrow (strings.size () > 0,
                 dealii::ExcMessage ("Fissile material IDs must be inserted for eigen problems"));
    // std::set<int> fissile_ids;
    for (int i=0; i<strings.size(); ++i)
      fissile_ids_.insert (std::atoi(strings[i].c_str())-1);
  }
  prm.leave_subsection ();

  for (int m=0; m<n_material_; ++m) {
    if (fissile_ids_.find(m)!=fissile_ids_.end())
      is_material_fissile_[m] = true;
    else
      is_material_fissile_[m] = false;
  }
  AssertThrow (!is_material_fissile_.empty (),
      dealii::ExcMessage ("Please specify at least one valid ID for fissile materials"));
  if (n_group_>1) {
    prm.enter_subsection ("chi, group=1 to G");
    {
      for (int m=0; m<n_material_;++m) {
        std::vector<double> tmp(n_group_);
        if (is_material_fissile_[m]) {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings =
              dealii::Utilities::split_string_list (prm.get(os.str()));
          AssertThrow (strings.size () == n_group_,
              dealii::ExcMessage ("n_group_ is not equal to group number of chi"));
          std::vector<double> tmp(n_group_);
          for (int g=0; g<n_group_; ++g)
            tmp[g] = std::atof (strings[g].c_str ());
        }
        chi_[m] = tmp;
      }
    }
    prm.leave_subsection ();

    prm.enter_subsection ("nu_sigf, group=1 to G");
    {
      for (int m=0; m<n_material_;++m) {
        std::vector<double> tmp(n_group_);
        if (is_material_fissile_[m]) {
          std::ostringstream os;
          os << "material " << m + 1;
          std::vector<std::string> strings =
              dealii::Utilities::split_string_list (prm.get (os.str ()));
          AssertThrow (strings.size () == n_group_,
              dealii::ExcMessage ("n_group_ is not equal to group number of nusigf"));
          std::vector<double> tmp(n_group_);
          for (int g=0; g<n_group_; ++g)
            tmp[g] = std::atof (strings[g].c_str ());
        }
        nu_sigf_[m] = tmp;
      }
    }
    prm.leave_subsection ();
  } else {
    // the following is for one-group
    prm.enter_subsection ("one-group chi");
    {
      std::vector<std::string> strings =
          dealii::Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material_,
          dealii::ExcMessage("One-group chi should have n_material_ entries"));
      std::vector<double> tmp_sigt (n_material_);
      for (int m=0; m<n_material_; ++m) {
        std::vector<double> tmp {is_material_fissile_[m] ?
            std::atof (strings[m].c_str()) : 0.0};
        chi_[m] = tmp;
      }
    }
    prm.leave_subsection ();
    prm.enter_subsection ("one-group nu_sigf");
    {
      std::vector<std::string> strings = dealii::Utilities::split_string_list (prm.get("values"));
      AssertThrow (strings.size()==n_material_,
          dealii::ExcMessage("One-group sigma_t should have n_material_ entries"));
      std::vector<double> tmp_sigt (n_material_);
      for (int m=0; m<n_material_; ++m)
      {
        std::vector<double> tmp {is_material_fissile_[m] ?
            std::atof (strings[m].c_str()) : 0.0};
        nu_sigf_[m] = tmp;
      }
    }
    prm.leave_subsection ();
  }
  for (int m=0; m<n_material_; ++m)
  {
    dealii::FullMatrix<double> tmp (n_group_, n_group_);
    if (is_material_fissile_[m])
      for (int gin=0; gin<n_group_; ++gin)
        for (int g=0; g<n_group_; ++g)
          tmp(gin,g) = chi_[m][g] * nu_sigf_[m][gin];

    fiss_transfer_[m] = tmp;
    fiss_transfer_per_ster_[m] = tmp;
    fiss_transfer_per_ster_[m] *= bconst::kInvFourPi;
  }
}

std::unordered_map<int, std::vector<double>>
Materials::GetSigT () const {
  return sigt_;
}

std::unordered_map<int, std::vector<double>>
Materials::GetInvSigT () const {
  return inv_sigt_;
}

std::unordered_map<int, std::vector<double>> Materials::GetQ() const {
  return q_;
}

std::unordered_map<int, std::vector<double>>
Materials::GetQPerSter() const {
  return q_per_ster_;
}

std::unordered_map<int, dealii::FullMatrix<double>>
Materials::GetSigS() const {
  return sigs_;
}

std::unordered_map<int, dealii::FullMatrix<double>>
Materials::GetSigSPerSter () const {
  return sigs_per_ster_;
}

std::unordered_map<int, dealii::FullMatrix<double>>
Materials::GetFissTransfer () const {
  return fiss_transfer_;
}

std::unordered_map<int, dealii::FullMatrix<double>>
Materials::GetFissTransferPerSter () const {
  return fiss_transfer_per_ster_;
}

std::unordered_map<int, std::vector<double>> Materials::GetNuSigf () const {
  return nu_sigf_;
}

std::unordered_map<int, bool> Materials::GetFissileIDMap () const {
  return is_material_fissile_;
}
