#include "../../src/aqdata/aq_base.h"
#include "../test_utilities.h"

void SetupParameters (dealii::ParameterHandler &prm) {
  prm.declare_entry ("have reflective BC", "false", dealii::Patterns::Bool(), "");
  prm.declare_entry ("transport model", "regular", dealii::Patterns::Selection("regular|ep"), "");
  prm.declare_entry ("angular quadrature order", "4", dealii::Patterns::Integer (), "");
  prm.declare_entry ("angular quadrature name", "gl",
                     dealii::Patterns::Selection ("gl"), "");
  prm.declare_entry ("number of groups", "1", dealii::Patterns::Integer (), "");
}

template<int dim>
void Test (dealii::ParameterHandler &prm) {
  std::unique_ptr<AQBase<dim>> gl_ptr =
  	  std::unique_ptr<AQBase<dim>> (new AQBase<dim>(prm));
  gl_ptr->MakeAQ ();
  auto wi = gl_ptr->GetAQWeights ();
  auto omega_i = gl_ptr->GetAQDirs ();
  for (int i=0; i<wi.size(); ++i)
  {
    dealii::deallog << "Weight: " << wi[i] << "; Omega: ";
    for (int j=0; j<dim; ++j)
      dealii::deallog << omega_i[i][j] << " ";
    dealii::deallog << std::endl;
  }
}

int main () {
  dealii::ParameterHandler prm;
  SetupParameters (prm);
  testing::InitLog ();
  dealii::deallog << "AQ for Gauss-Legendre S4 (1D)" << std::endl;
  // 1D test for non-ep equations
  dealii::deallog.push ("1D");
  Test<1> (prm);
  dealii::deallog.pop ();

  dealii::deallog << std::endl << std::endl << std::endl;

  prm.set ("transport model", "ep");
  dealii::deallog.push ("1D EP");
  Test<1> (prm);
  dealii::deallog.pop ();

  return 0;
}
