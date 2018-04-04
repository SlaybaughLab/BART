#include "../problem_definition.h"

#include "gtest/gtest.h"

class ProblemDefinitionTest : public ::testing::Test {
 protected:
  void SetUp ();

  void DeclareParamsTestLocalPrm ();

  void DeclareParamsTestGlobalPrm ();

  dealii::ParameterHandler local_prm_;
};

void ProblemDefinitionTest::SetUp () {
  bparams::DeclareParameters (local_prm_);
  // setup parameters without default values
  local_prm_.set ("reflective boundary names", "xmin");
  local_prm_.set ("x, y, z max values of boundary locations", "1.0, 2.0");
  local_prm_.set ("number of cells for x, y, z directions", "1, 3, 2");

  // setup global prm
  bparams::DeclareParameters ();
  bparams::GlobPrm.set ("reflective boundary names", "xmin");
  bparams::GlobPrm.set ("x, y, z max values of boundary locations", "1.0, 2.0, 1.");
  bparams::GlobPrm.set ("number of cells for x, y, z directions", "1, 3, 2");
}

void ProblemDefinitionTest::DeclareParamsTestLocalPrm () {
  EXPECT_EQ (local_prm_.get_integer("problem dimension"), 2);
  EXPECT_EQ (local_prm_.get("transport model"), "none");
  EXPECT_EQ (local_prm_.get("ho linear solver name"), "cg");
  EXPECT_EQ (local_prm_.get("ho preconditioner name"), "amg");
  EXPECT_EQ (local_prm_.get_double("ho ssor factor"), 1.0);
  EXPECT_EQ (local_prm_.get("nda linear solver name"), "none");
  EXPECT_EQ (local_prm_.get("nda preconditioner name"), "jacobi");
  EXPECT_EQ (local_prm_.get_double("nda ssor factor"), 1.0);
  EXPECT_EQ (local_prm_.get("angular quadrature name"), "none");
  EXPECT_EQ (local_prm_.get_integer("angular quadrature order"), 4);
  EXPECT_EQ (local_prm_.get_integer("number of groups"), 1);
  EXPECT_EQ (local_prm_.get_integer("thermal group boundary"), 0);
  EXPECT_EQ (local_prm_.get("ho spatial discretization"), "cfem");
  EXPECT_EQ (local_prm_.get("nda spatial discretization"), "cfem");
  EXPECT_EQ (local_prm_.get_bool("do eigenvalue calculations"), false);
  EXPECT_EQ (local_prm_.get_bool("do nda"), false);
  EXPECT_EQ (local_prm_.get_bool("have reflective BC"), false);
  EXPECT_EQ (local_prm_.get("reflective boundary names"), "xmin");
  EXPECT_EQ (local_prm_.get_integer("finite element polynomial degree"), 1);
  EXPECT_EQ (local_prm_.get_integer("uniform refinements"), 0);
  EXPECT_EQ (local_prm_.get("x, y, z max values of boundary locations"), "1.0, 2.0");
  EXPECT_EQ (local_prm_.get("number of cells for x, y, z directions"), "1, 3");
  EXPECT_EQ (local_prm_.get_integer("number of materials"), 1);
  EXPECT_EQ (local_prm_.get_bool("do print angular quadrature info"), true);
  EXPECT_EQ (local_prm_.get_bool("is mesh generated by deal.II"), true);
  EXPECT_EQ (local_prm_.get("output file name base"), "solu");
  EXPECT_EQ (local_prm_.get("mesh file name"), "mesh.msh");
}

void ProblemDefinitionTest::DeclareParamsTestGlobalPrm () {
  EXPECT_EQ (bparams::GlobPrm.get_integer("problem dimension"), 2);
  EXPECT_EQ (bparams::GlobPrm.get("transport model"), "none");
  EXPECT_EQ (bparams::GlobPrm.get("ho linear solver name"), "cg");
  EXPECT_EQ (bparams::GlobPrm.get("ho preconditioner name"), "amg");
  EXPECT_EQ (bparams::GlobPrm.get_double("ho ssor factor"), 1.0);
  EXPECT_EQ (bparams::GlobPrm.get("nda linear solver name"), "none");
  EXPECT_EQ (bparams::GlobPrm.get("nda preconditioner name"), "jacobi");
  EXPECT_EQ (bparams::GlobPrm.get_double("nda ssor factor"), 1.0);
  EXPECT_EQ (bparams::GlobPrm.get("angular quadrature name"), "none");
  EXPECT_EQ (bparams::GlobPrm.get_integer("angular quadrature order"), 4);
  EXPECT_EQ (bparams::GlobPrm.get_integer("number of groups"), 1);
  EXPECT_EQ (bparams::GlobPrm.get_integer("thermal group boundary"), 0);
  EXPECT_EQ (bparams::GlobPrm.get("ho spatial discretization"), "cfem");
  EXPECT_EQ (bparams::GlobPrm.get("nda spatial discretization"), "cfem");
  EXPECT_EQ (bparams::GlobPrm.get_bool("do eigenvalue calculations"), false);
  EXPECT_EQ (bparams::GlobPrm.get_bool("do nda"), false);
  EXPECT_EQ (bparams::GlobPrm.get_bool("have reflective BC"), false);
  EXPECT_EQ (bparams::GlobPrm.get("reflective boundary names"), "xmin");
  EXPECT_EQ (bparams::GlobPrm.get_integer("finite element polynomial degree"), 1);
  EXPECT_EQ (bparams::GlobPrm.get_integer("uniform refinements"), 0);
  EXPECT_EQ (bparams::GlobPrm.get("x, y, z max values of boundary locations"), "1.0, 2.0");
  EXPECT_EQ (bparams::GlobPrm.get("number of cells for x, y, z directions"), "1, 3");
  EXPECT_EQ (bparams::GlobPrm.get_integer("number of materials"), 1);
  EXPECT_EQ (bparams::GlobPrm.get_bool("do print angular quadrature info"), true);
  EXPECT_EQ (bparams::GlobPrm.get_bool("is mesh generated by deal.II"), true);
  EXPECT_EQ (bparams::GlobPrm.get("output file name base"), "solu");
  EXPECT_EQ (bparams::GlobPrm.get("mesh file name"), "mesh.msh");
}

TEST_F (ProblemDefinitionTest, DeclareParamsLocalTest) {
  DeclareParamsTestLocalPrm ();
}

TEST_F (ProblemDefinitionTest, DeclareParamsGlobalTest) {
  DeclareParamsTestGlobalPrm ();
}
