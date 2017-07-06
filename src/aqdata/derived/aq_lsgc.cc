#include <deal.II/base/quadrature_lib.h>

#include <iostream>

#include "../../../include/aqdata/derived/aq_lsgc.h"

template <int dim>
AQLSGC<dim>::AQLSGC (ParameterHandler &prm)
:
AQBase<dim> (prm)
{
}

template <int dim>
AQLSGC<dim>::~AQLSGC ()
{
}

template <int dim>
void AQLSGC<dim>::produce_angular_quad ()
{
  AssertThrow (dim>=2,
               ExcMessage("1D is not implemented"));
  AssertThrow (this->n_azi%2==0,
               ExcMessage("SN order must be even numbers"));
  QGauss<1> mu_quad (this->n_azi);

  this->total_angle =
  (dim==2 ? 8.0*this->pi :
   4.0*this->pi*(this->transport_model_name=="ep" ? 2.0 : 1.0));
  unsigned int n_total_azi = ((dim==3 && this->transport_model_name!="ep")?
                              this->n_azi : this->n_azi / 2);
  for (unsigned int i=0; i<n_total_azi; ++i)
  {
    double mu = mu_quad.point(i)[0] * 2.0 - 1.0;
    unsigned int n_level = ((i<this->n_azi/2?4*(i+1):4*(this->n_azi-i))/
                            ((dim==2&&this->transport_model_name=="ep")?2:1));
    double dphi = 2.0 * this->pi / (n_level*(this->transport_model_name=="ep"?2.0:1.0));
    double w_pt = mu_quad.weight(i) * this->total_angle / n_level;
    for (unsigned int j=0; j<n_level; ++j)
    {
      Tensor<1, dim> omega;
      double phi = (j + 0.5) * dphi;
      omega[0] = std::sqrt (1.0 - mu * mu) * cos (phi);
      omega[1] = std::sqrt (1.0 - mu * mu) * sin (phi);
      if (dim==3)
        omega[2] = mu;
      this->wi.push_back (w_pt);
      this->omega_i.push_back (omega);
    }
  }

  this->n_dir = ((dim==2) ? this->n_azi * (this->n_azi + 2) / (this->transport_model_name=="ep"?4:2) :
                 this->n_azi * (this->n_azi + 2) / (this->transport_model_name=="ep"?2:1));

  AssertThrow (this->n_dir==this->wi.size(),
               ExcMessage("calculated number of angles should be the same as number of angular weights"));
  AssertThrow (this->n_dir==this->omega_i.size(),
               ExcMessage("calculated number of angles should be the same as number of angles"));

  this->n_total_ho_vars = this->n_dir * this->n_group;
  // estimate tensor norm to do penalty method for EP
  if (this->transport_model_name=="ep" &&
      this->discretization=="dfem")
    for (unsigned int i=0; i<this->n_dir; ++i)
    {
      Tensor<2, dim> tensor_tmp = outer_product(this->omega_i[i], this->omega_i[i]);
      this->tensor_norms.push_back(tensor_tmp.norm());
    }
}

template class AQLSGC<2>;
template class AQLSGC<3>;
