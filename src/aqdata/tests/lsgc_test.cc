#include "../lsgc.h"

#include "gtest/gtest.h"
#include "../../test_helpers/bart_test_helper.h"

class LSGCTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    prm.declare_entry ("have reflective BC", "false",
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
  void AQDataTest();

  template<int dim>
  void InitRefBCAndOutput();
  
  dealii::ParameterHandler prm;
};

template <int dim>
void LSGCTest::AQDataTest() {
  std::unique_ptr<AQBase<dim>> lsgc_ptr =
      std::unique_ptr<AQBase<dim>> (new LSGC<dim>(prm));
  lsgc_ptr->MakeAQ();
  auto wi = lsgc_ptr->GetAQWeights();
  auto omega_i = lsgc_ptr->GetAQDirs();
  for (unsigned int i=0; i<wi.size(); ++i)
  {
    dealii::deallog << "Weight: " << wi[i] << "; Omega: ";
    for (int j=0; j<dim; ++j)
      dealii::deallog << omega_i[i][j] << " ";
    dealii::deallog << std::endl;
  }
}

template<int dim>
void LSGCTest::InitRefBCAndOutput () {
  // Initializes BC reflection directions and outputs them to the deallog
  std::unique_ptr<AQBase<dim>> aq_base_ptr =
      std::unique_ptr<AQBase<dim>> (new LSGC<dim>(prm));
  aq_base_ptr->MakeAQ();
  std::map<std::pair<int, int>, int > reflection_dirs =
      aq_base_ptr->GetRefDirInd();

  // Print values to the deallog:
  for (auto const& entry : reflection_dirs) {
    dealii::deallog << "Boundary ID: " << entry.first.first
                    << ";Current dir: " << entry.first.second
                    << ";Reflect dir: " << entry.second
                    << std::endl;
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
  std::string filename = "lsgc_2d_refl";
  prm.set("have reflective BC", "true");

  btest::GoldTestInit(filename);
  
  InitRefBCAndOutput<2>();
  
  btest::GoldTestRun(filename);
}
