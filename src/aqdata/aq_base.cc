/*
 Class header file
 */
#include "aq_base.h"

/*
 STL header files
 */
#include <fstream>
#include <iomanip>
/*
 deal.II header files
 */
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>

template <int dim>
AQBase<dim>::AQBase (dealii::ParameterHandler &prm)
    : k_pi(dealii::numbers::PI),
      have_reflective_bc_(prm.get_bool("have reflective BC")),
      transport_model_name_(prm.get("transport model")),
      n_azi_(prm.get_integer("angular quadrature order")),
      n_group_(prm.get_integer("number of groups")),
      aq_name_(prm.get("angular quadrature name")) {}

template <int dim>
AQBase<dim>::~AQBase () {}

template <int dim>
void AQBase<dim>::MakeAQ () {
  ProduceAQ ();
  InitCompInd ();
  InitRefBCInd ();
}

template <int dim>
void AQBase<dim>::InitRefBCInd () {
  // Note: here we assume square domain and assume user either
  // uses deal.II generated mesh or mesh from gmsh with proper
  // boundary IDs setup: {0,1,2,3,4,5} for {xmin,xmax,ymin,ymax,zmin,zmax}
  if (have_reflective_bc_) {
    Assert (dim>1,
            dealii::ExcNotImplemented());
    std::vector<dealii::Tensor<1, dim>> bnv;
    bnv.resize (2*dim);
    // All boundary normal vectors are assume to be parallel to axes
    // Then, only one component in each normal vector is nonzero
    bnv[0][0] = -1.0;
    bnv[1][0] = 1.0;
    if (dim>=2) {
      bnv[2][1] = -1.0;
      bnv[3][1] = 1.0;
      if (dim>=3) {
        bnv[4][2] = -1.0;
        bnv[5][2] = 1.0;
      }
    }
    for (int i=0; i<2*dim; ++i)
      for (int i_dir=0; i_dir<n_dir_; ++i_dir) {
        dealii::Tensor<1, dim> out_angle =
            (omega_i_[i_dir] - 2.0 * (bnv[i]*omega_i_[i_dir]) * bnv[i]);
        //(omega_i_[i_dir] *
        // (1.0 - 2.0 * (bnv[i] * omega_i_[i_dir])));
        for (int r_dir=0; r_dir<n_dir_; ++r_dir) {
          dealii::Tensor<1, dim> d_dir = out_angle;
          dealii::Tensor<1, dim> d_minus_dir = out_angle;
          d_dir -= omega_i_[r_dir];
          d_minus_dir += omega_i_[r_dir];
          // Use caution about even parity using this std::map
          if ((transport_model_name_=="ep" &&
               (d_dir.norm ()<1.0e-13 || d_minus_dir.norm ()<1.0e-13)) ||
              (transport_model_name_!="ep" && d_dir.norm ()<1.0e-13))
            reflective_direction_index_[std::make_pair (i, i_dir)] = r_dir;
        }
      }
  }
}

template <int dim>
void AQBase<dim>::ProduceAQ () {
  // Only provide 1D by default
  AssertThrow (dim==1,
               dealii::ExcMessage("Only 1D is implemented in AQBase<dim>"));
  n_dir_ = n_azi_ / (transport_model_name_=="ep" ? 2 : 1);
  total_angle_ = transport_model_name_=="ep" ? 4.0 : 2.0;
  dealii::QGauss<1> mu_quad (n_azi_);
  for (int i=0; i<n_dir_; ++i) {
    dealii::Tensor <1,dim> omega;
    omega[0] = mu_quad.point(i)[0]*2.0-1.0;
    omega_i_.push_back (omega);
    wi_.push_back (mu_quad.weight(i)*total_angle_);
  }
  n_total_ho_vars_ = n_dir_ * n_group_;
}

template <int dim>
void AQBase<dim>::InitCompInd () {
  // initialize the map from (group, direction) to component indices
  int ind = 0;
  for (int g=0; g<n_group_; ++g)
    for (int i_dir=0; i_dir<n_dir_; ++i_dir)
    {
      std::pair<int, int> key (g, i_dir);
      component_index_[key] = ind;
      inverse_component_index_[ind] = key;
      ind += 1;
    }
}

template <int dim>
void AQBase<dim>::PrintAQ () {
  std::ofstream quadr;
  quadr.open("aq.txt");
  quadr << "transport model: " << transport_model_name_
        << "; quadrature name: " << ProduceAQName () << std::endl;
  quadr << "Dim = " << dim << ", SN order = " << n_azi_ << std::endl;
  quadr << "Weights | Omega_x | Omega_y | mu" << std::endl;
  for (int i=0; i<omega_i_.size(); ++i) {
    quadr << std::fixed << std::setprecision (15)
          << wi_[i] << "  " << omega_i_[i][0] << "  ";
    if (dim>1) {
      quadr << omega_i_[i][1] << "  ";
      double mu = dim>2 ? omega_i_[i][2] :
          std::sqrt(1.-(std::pow(omega_i_[i][0],2.)
                        +std::pow(omega_i_[i][1],2)));
      quadr << mu;
    }
    quadr << std::endl;
  }
  quadr.close ();
}

template <int dim>
std::string AQBase<dim>::ProduceAQName () {
  AssertThrow (aq_name_.size()>0,
               dealii::ExcMessage("aq name has to be assigned"));
  // ToDo: more quadrature name producers
  if (aq_name_=="lsgc")
    return "Level Symmetric Gauss Chebyshev";
}

//public member functions to retrieve private and protected variables
template <int dim>
std::map<std::pair<int, int>, int> AQBase<dim>::GetCompInd () {
  return component_index_;
}

template <int dim>
std::unordered_map<int, std::pair<int, int>> AQBase<dim>::GetInvCompInd () {
  return inverse_component_index_;
}

template <int dim>
std::map<std::pair<int, int>, int> AQBase<dim>::GetRefDirInd () {
  return reflective_direction_index_;
}

template <int dim>
int AQBase<dim>::GetSnOrder () {
  return n_azi_;
}

template <int dim>
int AQBase<dim>::GetNDir () {
  return n_dir_;
}

template <int dim>
int AQBase<dim>::GetNTotalHOVars () {
  return n_total_ho_vars_;
}

template <int dim>
std::vector<double> AQBase<dim>::GetAQWeights () {
  return wi_;
}

template <int dim>
std::vector<dealii::Tensor<1, dim>> AQBase<dim>::GetAQDirs () {
  return omega_i_;
}

template class AQBase<1>;
template class AQBase<2>;
template class AQBase<3>;
