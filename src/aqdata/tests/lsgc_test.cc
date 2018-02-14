#include "gtest/gtest.h"

#include "../lsgc.h"
#include "../../test_helpers/bart_test_helper.h"

class LSGCTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    prm.declare_entry ("have reflective BC", "false",
                       dealii::Patterns::Bool(), "");
    prm.declare_entry ("transport model", "ep",
                       dealii::Patterns::Selection("ep"), "");
    prm.declare_entry ("angular quadrature order", "4",
                       dealii::Patterns::Integer (), "");
    prm.declare_entry ("angular quadrature name", "lsgc",
                       dealii::Patterns::Selection ("lsgc"), "");
    prm.declare_entry ("number of groups", "1",
                       dealii::Patterns::Integer (), "");
  }
  template <int dim>
  void AQDataTest() {
    std::unique_ptr<AQBase<dim>> lsgc_ptr =
        std::unique_ptr<AQBase<dim>> (new LSGC<2>(prm));
    lsgc_ptr->make_aq ();
    auto wi = lsgc_ptr->get_angular_weights ();
    auto omega_i = lsgc_ptr->get_all_directions ();
    for (unsigned int i=0; i<wi.size(); ++i)
    {
      dealii::deallog << "Weight: " << wi[i] << "; Omega: ";
      for (int j=0; j<dim; ++j)
        dealii::deallog << omega_i[i][j] << " ";
      dealii::deallog << std::endl;
    }
  }
  
  dealii::ParameterHandler prm;
};

TEST_F(LSGCTest, 2DLSGCTest) {
  std::string filename = "2dlsgc";
  btest::GlobalBartTestHelper().OpenLog(filename);
  
  AQDataTest<2>();
  
  btest::GlobalBartTestHelper().CloseLog();
  ASSERT_TRUE(btest::GlobalBartTestHelper().GoldTest(filename));  
}
