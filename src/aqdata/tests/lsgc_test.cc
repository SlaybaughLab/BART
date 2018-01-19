#include "gtest/gtest.h"
#include "../lsgc.h"

#include <deal.II/base/logstream.h>

class LSGCTest : public ::testing::Test
{
 protected:
  virtual void SetUp() {
    dealii::ParameterHandler prm;
    setup_parameters(prm);
    LSGC<2> lsgc_2_(prm);
    lsgc_2_.make_aq();
  }

  void setup_parameters (dealii::ParameterHandler &prm)
  {
    prm.declare_entry ("have reflective BC", "false", dealii::Patterns::Bool(), "");
    prm.declare_entry ("transport model", "ep", dealii::Patterns::Selection("ep"), "");
    prm.declare_entry ("angular quadrature order", "4", dealii::Patterns::Integer (), "");
    prm.declare_entry ("angular quadrature name", "lsgc",
                       dealii::Patterns::Selection ("lsgc"), "");
    prm.declare_entry ("number of groups", "1", dealii::Patterns::Integer (), "");
  }

  LSGC<2> *lsgc_2_;
 
}; 
  

TEST (LSGCTest, 2DTest) {
  ASSERT_EQ(0, 0);
}
