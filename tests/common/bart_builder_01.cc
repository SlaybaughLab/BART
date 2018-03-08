#include "../test_utilities.h"
#include "../../src/common/bart_builder.h"

void DeclareParams (dealii::ParameterHandler &prm) {
  prm.declare_entry ("do nda", "true", dealii::Patterns::Bool(), "");
  prm.declare_entry ("finite element polynomial degree", "1",
                     dealii::Patterns::Integer(), "");
  prm.declare_entry ("ho spatial discretization", "",
                     dealii::Patterns::Anything(), "");
  prm.declare_entry ("nda spatial discretization", "",
                     dealii::Patterns::Anything(), "");
}

template <int dim>
void Test (dealii::ParameterHandler &prm) {
  dealii::deallog.push (dealii::Utilities::int_to_string(dim)+"D");
  BARTBuilder<dim> builders (prm);

  prm.set ("ho spatial discretization", "cfem");
  prm.set ("nda spatial discretization", "cfem");
  builders.SetParams (prm);
  std::vector<dealii::FiniteElement<dim, dim>*> fe_ptrs;
  builders.BuildFESpaces (fe_ptrs);
  dealii::deallog << "HO FE is " << fe_ptrs.front()->get_name() << std::endl;
  dealii::deallog << "NDA FE is " << fe_ptrs.back()->get_name() << std::endl;

  prm.set ("ho spatial discretization", "dfem");
  prm.set ("nda spatial discretization", "dfem");
  builders.SetParams (prm);
  fe_ptrs.clear ();
  builders.BuildFESpaces (fe_ptrs);
  dealii::deallog << "HO FE is " << fe_ptrs.front()->get_name() << std::endl;
  dealii::deallog << "NDA FE is " << fe_ptrs.back()->get_name() << std::endl;

  prm.set ("ho spatial discretization", "dfem");
  prm.set ("nda spatial discretization", "cmfd");
  builders.SetParams (prm);
  fe_ptrs.clear ();
  builders.BuildFESpaces (fe_ptrs);
  dealii::deallog << "HO FE is " << fe_ptrs.front()->get_name() << std::endl;
  dealii::deallog << "NDA FE is " << fe_ptrs.back()->get_name() << std::endl;

  prm.set ("ho spatial discretization", "dfem");
  prm.set ("nda spatial discretization", "rtk");
  builders.SetParams (prm);
  fe_ptrs.clear ();
  builders.BuildFESpaces (fe_ptrs);
  dealii::deallog << "HO FE is " << fe_ptrs.front()->get_name() << std::endl;
  dealii::deallog << "NDA FE is " << fe_ptrs.back()->get_name() << std::endl;

  dealii::deallog.pop ();
}

int main () {
  testing::InitLog ();
  dealii::ParameterHandler prm;
  DeclareParams (prm);
  Test<2> (prm);
  return 0;
}
