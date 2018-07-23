#include "../lsgc.h"

#include <sstream>

#include "gtest/gtest.h"
#include "../../test_helpers/bart_test_helper.h"

class LSGCTest : public ::testing::Test {
 protected:
  void SetUp() override;

  template <int dim>
  void AQDataTest();

  template<int dim>
  void InitRefBCAndOutput();

  template<int dim>
  void InitRefBCAndCheck();
  
  dealii::ParameterHandler prm;
};

void LSGCTest::SetUp() {
    prm.declare_entry ("have reflective boundary", "false",
                       dealii::Patterns::Bool(), "");
    prm.declare_entry ("angular quadrature order", "4",
                       dealii::Patterns::Integer (), "");
    prm.declare_entry ("angular quadrature name", "lsgc",
                       dealii::Patterns::Selection ("lsgc"), "");
    prm.declare_entry ("number of groups", "1",
                       dealii::Patterns::Integer (), "");
    prm.declare_entry ("transport model", "regular",
                       dealii::Patterns::Selection("regular|ep"), "");
}

template <int dim>
void LSGCTest::AQDataTest() {
  std::unique_ptr<LSGC<dim>> lsgc_ptr = std::make_unique<LSGC<dim>>(prm);
  lsgc_ptr->MakeAQ();
  auto wi = lsgc_ptr->GetAQWeights();
  auto omega_i = lsgc_ptr->GetAQDirs();
  for (int i=0; i<wi.size(); ++i)
  {
    dealii::deallog << "Weight: " << wi[i] << "; Omega: ";
    for (int j=0; j<dim; ++j)
      dealii::deallog << omega_i[i][j] << " ";
    dealii::deallog << std::endl;
  }
}

template<int dim>
void LSGCTest::InitRefBCAndCheck() {
  std::unique_ptr<LSGC<dim>> aq_base_ptr = std::make_unique<LSGC<dim>>(prm);
  aq_base_ptr->MakeAQ();
  // Get omegas and reflection map
  std::vector<dealii::Tensor<1, dim>> omegas = aq_base_ptr->GetAQDirs();
  std::map<std::pair<int, int>, int > reflection_map =
      aq_base_ptr->GetRefDirInd();

  for (auto const& mapping : reflection_map) {
    // Unroll this data structure
    int boundary = mapping.first.first;
    int direction = mapping.first.second;
    int reflection = mapping.second;

    switch (boundary / 2) {
      case 0: {
        //x-direction
        EXPECT_FLOAT_EQ(omegas[direction][0], -omegas[reflection][0]);
        break;
      }
      case 1: {
        //y-direction
        EXPECT_FLOAT_EQ(omegas[direction][1], -omegas[reflection][1]);
        break;
      }
      case 2: {
        //z-direction
        EXPECT_FLOAT_EQ(omegas[direction][2], -omegas[reflection][2]);
        break;
      }
      default: {
        assert(false);
      }
    }
  } 
}

TEST_F (LSGCTest, LSGC_2d_EpTest) {
  std::string filename = "lsgc_ep_2d";
  btest::GoldTestInit(filename);
  prm.set ("transport model", "ep"); 
  AQDataTest<2>();
  btest::GoldTestRun(filename);
}

TEST_F (LSGCTest, LSGC_3d_EpTest) {
  std::string filename = "lsgc_ep_3d";
  btest::GoldTestInit(filename);
  prm.set ("transport model", "ep");
  AQDataTest<3>(); 
  btest::GoldTestRun(filename);
}

TEST_F (LSGCTest, LSGC_2d_Test) {
  std::string filename = "lsgc_2d";
  btest::GoldTestInit(filename);
  AQDataTest<2>();
  btest::GoldTestRun(filename);
}

TEST_F (LSGCTest, LSGC_3d_Test) {
  std::string filename = "lsgc_3d";
  btest::GoldTestInit(filename);
  AQDataTest<3>();
  btest::GoldTestRun(filename);
}

TEST_F (LSGCTest, AQBase2DReflDir) {
  prm.set("have reflective boundary", "true");
  InitRefBCAndCheck<2>();
}

TEST_F (LSGCTest, AQBase3DReflDir) {
  prm.set("have reflective boundary", "true");
  InitRefBCAndCheck<3>();
}

TEST_F (LSGCTest, LSGC_2d_PrintAQ) {
  LSGC<2> test_lsgc(prm);
  test_lsgc.MakeAQ();
  
  std::ostringstream output_string_stream;
  test_lsgc.PrintAQ(&output_string_stream);
  
  std::string output = "transport model: regular; output_streamature name: "
                       "Level Symmetric Gauss Chebyshev\n"
                       "Dim = 2, SN order = 4\n"
                       "Weights | Omega_x | Omega_y | mu\n"
                       "1.092818225999440  0.359474792477992  0.359474792477992  0.861136311594053\n"
                       "1.092818225999440  -0.359474792477992  0.359474792477992  0.861136311594053\n"
                       "1.092818225999440  -0.359474792477992  -0.359474792477992  0.861136311594053\n"
                       "1.092818225999440  0.359474792477992  -0.359474792477992  0.861136311594053\n"
                       "1.024387213795176  0.868846143426105  0.359887856222652  0.339981043584856\n"
                       "1.024387213795176  0.359887856222652  0.868846143426105  0.339981043584856\n"
                       "1.024387213795176  -0.359887856222652  0.868846143426105  0.339981043584856\n"
                       "1.024387213795176  -0.868846143426105  0.359887856222652  0.339981043584856\n"
                       "1.024387213795176  -0.868846143426105  -0.359887856222652  0.339981043584856\n"
                       "1.024387213795176  -0.359887856222653  -0.868846143426105  0.339981043584856\n"
                       "1.024387213795176  0.359887856222652  -0.868846143426105  0.339981043584856\n"
                       "1.024387213795176  0.868846143426105  -0.359887856222653  0.339981043584856\n";

  EXPECT_EQ(output, output_string_stream.str());
}


