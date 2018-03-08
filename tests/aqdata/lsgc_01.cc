#include "../../src/aqdata/lsgc.h"
#include "../test_utilities.h"

template <int dim>
void Test (dealii::ParameterHandler &prm) {
  std::unique_ptr<AQBase<dim>> lsgc_ptr =
  	  std::unique_ptr<AQBase<dim>> (new LSGC<dim>(prm));
  lsgc_ptr->MakeAQ ();
  auto wi = lsgc_ptr->GetAQWeights ();
  auto omega_i = lsgc_ptr->GetAQDirs ();
  for (int i=0; i<wi.size(); ++i) {
    dealii::deallog << "Weight: " << wi[i] << "; Omega: ";
    for (int j=0; j<dim; ++j)
      dealii::deallog << omega_i[i][j] << " ";
    dealii::deallog << std::endl;
  }
}

void SetupParameters (dealii::ParameterHandler &prm) {
  prm.declare_entry ("have reflective BC", "false", dealii::Patterns::Bool(), "");
  prm.declare_entry ("transport model", "ep", dealii::Patterns::Selection("ep"), "");
  prm.declare_entry ("angular quadrature order", "4", dealii::Patterns::Integer (), "");
  prm.declare_entry ("angular quadrature name", "lsgc",
                     dealii::Patterns::Selection ("lsgc"), "");
  prm.declare_entry ("number of groups", "1", dealii::Patterns::Integer (), "");
}

int main () {
  dealii::ParameterHandler prm;
  SetupParameters (prm);

  testing::InitLog ();

  dealii::deallog << "AQ for Even-Parity S4" << std::endl;
  // 2D test
  dealii::deallog.push ("2D");
  Test<2> (prm);
  dealii::deallog.pop ();

  dealii::deallog << std::endl << std::endl << std::endl;

  // 3D test
  dealii::deallog.push ("3D");
  Test<3> (prm);
  dealii::deallog.pop ();

  return 0;
}
