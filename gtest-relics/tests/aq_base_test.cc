#include "../aq_base.h"

// include gtest for basic testing functions
#include "gtest/gtest.h"

// include gmock for testing functions for STL containers
#include "gmock/gmock.h"

using ::testing:::ElementsAre;

template <int dim>
class AQDerivedMock : public AQBase<dim>
{
public:
  AQDerivedMock (dealii::ParameterHandler &prm);
  ~AQDerivedMock ();

  void produce_angular_quad ();
};

template <int dim>
AQDerivedMock<dim>::AQDerivedMock (dealii::ParameterHandler &prm)
    : AQBase<dim> (prm)
{}

template <int dim>
AQDerivedMock<dim>::~AQDerivedMock ()
{}

template <int dim>
void AQDerivedMock<dim>::produce_angular_quad ()
{
  // we use S2 quad for testing

  this->wi = (dim==2? std::vector<double> (4, this->k_pi) :
              std::vector<double> (8, this->k_pi));
  double omega = 0.5773502691896257;
  double omegas[4][2] = {{omega, omega},
                         {-omega, omega},
                         {-omega, -omega},
                         {omega, -omega}};
  this->omega_i = {dealii::Tensor<1,dim> (omegas[0]),
                   dealii::Tensor<1,dim> (omegas[1]),
                   dealii::Tensor<1,dim> (omegas[2]),
                   dealii::Tensor<1,dim> (omegas[3])};
  this->n_dir = this->omega_i.size();
  this->n_total_ho_vars = this->n_dir * this->n_group;
}

TEST (AQBaseTest, AQBaseTest2D)
{
  // set testing parameters for mock
  dealii::ParameterHandler prm;
  
  // mock relevant input strings
  std::string ref_bc = "true";
  std::string transport_model = "saaf";
  std::string num_grp = "2";
  std::string n_azi_angle = "2";
  std::string aq_name = "mock";
  
  // set parameter handler entries with mocking values
  prm.set ("have reflective BC", ref_bc);
  prm.set ("transport model", transport_model);
  prm.set ("angular quadrature order", n_azi_angle);
  prm.set ("number of groups", num_grp);
  prm.set ("angular quadrature name", aq_name);
  
  // TODO: make sure the following specifying dimensions is the
  // correct way
  std::unique_ptr<AQBase<2>> aq_mock =
      std::make_unique<AQBase<2>>(new AQDerivedMock(prm));
  aq_mock->make_aq ();
  
  const double k_pi = 3.14159265358979323846;
  
  // real and mock angular weights
  std::vector<double> real_wi (4, k_pi);
  std::vector<double> mock_wi = aq_mock->get_angular_weights ();

  // real and mock direction tensors
  double omega = 0.5773502691896257;
  double omegas[4][2] = {{omega, omega},
                         {-omega, omega},
                         {-omega, -omega},
                         {omega, -omega}};
  std::vector<dealii::Tensor<1, 2>> real_omega_i = {dealii::Tensor<1,dim> (omegas[0]),
                   							   	    dealii::Tensor<1,dim> (omegas[1]),
                   							        dealii::Tensor<1,dim> (omegas[2]),
                   							   	    dealii::Tensor<1,dim> (omegas[3])};
  std::vector<dealii::Tensor<1, 2>> mock_omega_i = aq_mock->get_all_directions ();
  
  // component mapping: (group, direction)->component
  // inverse mapping: component->(group, direction)
  std::map<std::pair<int, int>, int> real_cmp_idx;
  std::unordered_map<int, std::pair<int, int>> real_inv_cmp_idx;
  int idx = 0;
  for (int g=0; g<2; ++g)
  	for (int i_dir=0; i_dir<4; ++i_dir)
  	{
  	  real_cmp_idx[std::make_pair(g, i_dir)] = idx;
  	  real_inv_cmp_idx[idx++] = std::make_pair (g, i_dir);
  	}
  auto mock_cmp_idx = aq_mock->get_component_index ();
  auto mock_inv_cmp_idx = aq->get_inv_component_map ();
  // getting useful info from mock aq
  // use expect so all the failures will be reported
  EXPECT_EQ (2, aq_mock->get_sn_order());
  EXPECT_EQ (4, aq_mock->get_n_dir());
  EXPECT_EQ (8, aq_mock->get_n_total_vars());

  EXPECT_THAT (real_wi, ContainerEq(mock_wi));
  EXPECT_THAT (real_omega_i, ContainerEq(mock_omega_i));
  EXPECT_THAT (real_cmp_idx, ContainerEq(mock_cmp_idx));
  EXPECT_THAT (real_inv_cmp_idx, ContainerEq(mock_inv_cmp_idx));
}
