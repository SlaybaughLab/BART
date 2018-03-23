#include "../bart_builder.h"
#include "../problem_definition.h"

#include "gtest/gtest.h"
#include "../../test_helpers/bart_test_helper.h"

class BARTBuilderTest : public ::testing::Test {
 protected:
  void SetUp ();

  // FE builder test
  template <int dim>
  void FEBuilderTest ();

  // AQ builder test
  template <int dim>
  void AQBuilderTest ();

  dealii::ParameterHandler prm;
};

void BARTBuilderTest::SetUp () {
  bparams::DeclareParameters(prm);
  prm.set("do nda", "true");
}

template <int dim>
void BARTBuilderTest::FEBuilderTest () {
  BARTBuilder<dim> builders (prm);

  // set values to parameters
  prm.set ("ho spatial discretization", "cfem");
  prm.set ("nda spatial discretization", "cfem");
  builders.SetParams (prm);
  std::vector<dealii::FiniteElement<dim, dim>*> fe_ptrs;
  builders.BuildFESpaces (fe_ptrs);

  // testing for FE names
  EXPECT_EQ (fe_ptrs.front()->get_name(),
      "FE_Q<"+dealii::Utilities::int_to_string(dim)+">(1)");
  EXPECT_EQ (fe_ptrs.back()->get_name(),
      "FE_Q<"+dealii::Utilities::int_to_string(dim)+">(1)");

  // changing FE types
  prm.set ("ho spatial discretization", "dfem");
  prm.set ("nda spatial discretization", "dfem");
  builders.SetParams (prm);
  fe_ptrs.clear ();
  builders.BuildFESpaces (fe_ptrs);
  EXPECT_EQ (fe_ptrs.front()->get_name(),
      "FE_DGQ<"+dealii::Utilities::int_to_string(dim)+">(1)");
  EXPECT_EQ (fe_ptrs.back()->get_name(),
      "FE_DGQ<"+dealii::Utilities::int_to_string(dim)+">(1)");

  // changing NDA FE type
  prm.set ("nda spatial discretization", "cmfd");
  builders.SetParams (prm);
  fe_ptrs.clear ();
  builders.BuildFESpaces (fe_ptrs);
  EXPECT_EQ (fe_ptrs.back()->get_name(),
      "FE_DGQ<"+dealii::Utilities::int_to_string(dim)+">(0)");

  // changing NDA FE type
  prm.set ("nda spatial discretization", "rtk");
  builders.SetParams (prm);
  fe_ptrs.clear ();
  builders.BuildFESpaces (fe_ptrs);
  EXPECT_EQ (fe_ptrs.back()->get_name(),
      "FE_RaviartThomas<"+dealii::Utilities::int_to_string(dim)+">(1)");
}

TEST_F (BARTBuilderTest, FEBuilder2DTest) {
  FEBuilderTest<2> ();
}

TEST_F (BARTBuilderTest, FEBuilder3DTest) {
  FEBuilderTest<3> ();
}

template <int dim>
void BARTBuilderTest::AQBuilderTest () {
  BARTBuilder<dim> builders (prm);
  std::unique_ptr<AQBase<dim>> aq_ptr;
  std::string aq_name = (dim==1) ? "gl" : "lsgc";
  prm.set ("angular quadrature name", aq_name);
  prm.set ("angular quadrature order", "4");

  // call builder function to build aq
  builders.BuildAQ (prm, aq_ptr);
  aq_ptr->MakeAQ();
  auto wi = aq_ptr->GetAQWeights();
  auto omega_i = aq_ptr->GetAQDirs();

  // check output
  std::string filename = "aq_builder_"+dealii::Utilities::int_to_string(dim)+"d";
  btest::GoldTestInit(filename);
  for (unsigned int i=0; i<wi.size(); ++i)
  {
    dealii::deallog << "Weight: " << wi[i] << "; Omega: ";
    for (int j=0; j<dim; ++j)
      dealii::deallog << omega_i[i][j] << " ";
    dealii::deallog << std::endl;
  }
  btest::GoldTestRun(filename);
}

TEST_F (BARTBuilderTest, AQBuilder1DTest) {
  AQBuilderTest<1> ();
}

TEST_F (BARTBuilderTest, AQBuilder2DTest) {
  AQBuilderTest<2> ();
}

TEST_F (BARTBuilderTest, AQBuilder3DTest) {
  AQBuilderTest<3> ();
}
