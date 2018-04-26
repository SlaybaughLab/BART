#include "lsgc.h"

#include <deal.II/base/quadrature_lib.h>

template <int dim>
LSGC<dim>::LSGC (const dealii::ParameterHandler &prm)
:
AQBase<dim> (prm) {}

template <int dim>
LSGC<dim>::~LSGC () {}

template <int dim>
void LSGC<dim>::ProduceAQ () {
  AssertThrow (this->n_azi_%2==0,
               dealii::ExcMessage("SN order must be even numbers"));
  AssertThrow (dim>1,
               dealii::ExcMessage("LSGC only exists in multi-D"));

  dealii::QGauss<1> mu_quad (this->n_azi_);

  switch (dim) {
    case 2: {
      this->n_dir_ = this->n_azi_ * (this->n_azi_ + 2) /
          (this->transport_model_name_=="ep" ? 4:2);
      this->total_angle_ = 8.0 * this->k_pi;
      break;
    }
    case 3: {
      this->n_dir_ = this->n_azi_ * (this->n_azi_ + 2) /
          (this->transport_model_name_=="ep" ? 2:1);
      this->total_angle_ = 4.0 * this->k_pi *
          (this->transport_model_name_=="ep" ? 2.0 : 1.0);
      break;
    }
    default:
      break;
  }

  unsigned int n_total_azi = ((dim==3 && this->transport_model_name_!="ep")?
                              this->n_azi_ : this->n_azi_ / 2);
  for (unsigned int i=0; i<n_total_azi; ++i) {
    double mu = mu_quad.point(i)[0] * 2.0 - 1.0;
    unsigned int n_level = ((i<this->n_azi_/2?4*(i+1):4*(this->n_azi_-i))/
                            ((dim==2&&this->transport_model_name_=="ep")?2:1));
    double dphi = 2.0 * this->k_pi /
        (n_level*(this->transport_model_name_=="ep"?2.0:1.0));
    double w_pt = mu_quad.weight(i) * this->total_angle_ / n_level;
    for (unsigned int j=0; j<n_level; ++j) {
      dealii::Tensor<1, dim> omega;
      double phi = (j + 0.5) * dphi;
      omega[0] = std::sqrt (1.0 - mu * mu) * cos (phi);
      omega[1] = std::sqrt (1.0 - mu * mu) * sin (phi);
      if (dim==3)
        omega[2] = mu;
      this->wi_.push_back (w_pt);
      this->omega_i_.push_back (omega);
    }
  }

  AssertThrow (this->n_dir_==this->wi_.size(),
               dealii::ExcMessage("calculated number of angles should be the same as number of angular weights"));
  AssertThrow (this->n_dir_==this->omega_i_.size(),
               dealii::ExcMessage("calculated number of angles should be the same as number of angles"));

  this->n_total_ho_vars_ = this->n_dir_ * this->n_group_;
}

template class LSGC<2>;
template class LSGC<3>;
