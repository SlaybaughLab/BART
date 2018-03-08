#include "../../src/aqdata/aq_base.h"
#include "../test_utilities.h"


template <int dim>
class AQDerivedMock : public AQBase<dim> {
public:
  AQDerivedMock (dealii::ParameterHandler &prm);
  ~AQDerivedMock ();

  void ProduceAQ ();
  void InitCompInd ();
};

template <int dim>
AQDerivedMock<dim>::AQDerivedMock (dealii::ParameterHandler &prm)
    : AQBase<dim> (prm) {}

template <int dim>
AQDerivedMock<dim>::~AQDerivedMock () {}

template <int dim>
void AQDerivedMock<dim>::ProduceAQ () {
  dealii::deallog << "Mocking producing angular quadrature" << std::endl;

  this->wi_ = (dim==2 ? std::vector<double> (4, this->k_pi) :
               std::vector<double> (8, this->k_pi));
  this->omega_i_.resize (this->wi_.size());
  for (int i=0; i<this->wi_.size(); ++i)
    for (int j=0; j<dim; ++j)
      this->omega_i_[i][j] = i;
  this->n_dir_ = this->omega_i_.size();
  this->n_total_ho_vars_ = this->n_dir_ * this->n_group_;
}

template <int dim>
void AQDerivedMock<dim>::InitCompInd () {
  dealii::deallog << "Mocking initializing component index" << std::endl;
}

template <int dim>
void Test () {
  dealii::ParameterHandler prm;

  prm.declare_entry ("have reflective BC", "false", dealii::Patterns::Bool (), "");
  prm.declare_entry ("transport model", "mock", dealii::Patterns::Anything (), "");
  prm.declare_entry ("angular quadrature order", "2", dealii::Patterns::Integer (), "");
  prm.declare_entry ("angular quadrature name", "mock", dealii::Patterns::Anything (), "");
  prm.declare_entry ("number of groups", "2", dealii::Patterns::Integer (), "");

  std::unique_ptr<AQBase<dim>> aq_mock =
      std::unique_ptr<AQBase<dim>>(new AQDerivedMock<dim>(prm));
  aq_mock->MakeAQ ();
  auto omega_i = aq_mock->GetAQDirs ();
  auto wi = aq_mock->GetAQWeights ();
  dealii::deallog << "SN order: " << aq_mock->GetSnOrder () << "; ";
  dealii::deallog << "Total components: " << aq_mock->GetNTotalHOVars() << std::endl;
  // the first part is get_xxx functionality
  for (int i=0; i<wi.size(); ++i) {
    dealii::deallog << "Weight: " << wi[i] << "; Omega: ";
    for (int j=0; j<dim; ++j)
      dealii::deallog << omega_i[i][j] << " ";
    dealii::deallog << std::endl;
  }
}

int main () {
  testing::InitLog ();

  // put the test code here
  dealii::deallog.push ("2D");
  Test<2> ();
  dealii::deallog.pop ();
  dealii::deallog << std::endl
                  << "++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl << std::endl;
  dealii::deallog.push ("3D");
  Test<3> ();
  dealii::deallog.pop ();

  return 0;
}
