#include "../../src/aqdata/lsgc.h"

#include <fstream>

#include <deal.II/base/logstream.h>

template <int dim>
void test (dealii::ParameterHandler &prm)
{
  std::unique_ptr<AQBase<dim>> lsgc_ptr =
  	  std::unique_ptr<AQBase<dim>> (new LSGC<dim>(prm));
  lsgc_ptr->make_aq ();
  auto wi = lsgc_ptr->get_angular_weights ();
  auto omega_i = lsgc_ptr->get_all_directions ();
  for (int i=0; i<wi.size(); ++i)
  {
    dealii::deallog << "Weight: " << wi[i] << "; Omega: ";
    for (int j=0; j<dim; ++j)
      dealii::deallog << omega_i[i][j] << " ";
    dealii::deallog << std::endl;
  }
}

void setup_parameters (dealii::ParameterHandler &prm)
{
  prm.declare_entry ("have reflective BC", "false", dealii::Patterns::Bool(), "");
  prm.declare_entry ("transport model", "regular", dealii::Patterns::Selection("regular"), "");
  prm.declare_entry ("angular quadrature order", "4", dealii::Patterns::Integer (), "");
  prm.declare_entry ("angular quadrature name", "lsgc",
                     dealii::Patterns::Selection ("lsgc"), "");
  prm.declare_entry ("number of groups", "1", dealii::Patterns::Integer (), "");
}

int main ()
{
  dealii::ParameterHandler prm;
  setup_parameters (prm);
  std::string deallogname = "output";
  std::ofstream deallogfile (deallogname.c_str());
  dealii::deallog.attach (deallogfile, false);
  dealii::deallog << "AQ for Regular S4" << std::endl;
  // 2D test
  dealii::deallog.push ("2D");
  test<2> (prm);
  dealii::deallog.pop ();

  dealii::deallog << std::endl << std::endl << std::endl;

  // 3D test
  dealii::deallog.push ("3D");
  test<3> (prm);
  dealii::deallog.pop ();

  return 0;
}
