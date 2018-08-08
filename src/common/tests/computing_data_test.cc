#include "../computing_data.h"
#include "../problem_definition.h"
#include "../../test_helpers/bart_test_helper.h"

class ComputingDataTest
    :
    public btest::BARTParallelEnvironment {
 public:
  void SetUp ();

  dealii::ParameterHandler prm_;
};

void ComputingDataTest::SetUp() {
  // Initilize MPI
  this->MPIInit();
  // declare all entries for parameter handler
  bparams::DeclareParameters(prm_);

  prm_.set ("reflective boundary names", "xmin");
  prm_.set ("x, y, z max values of boundary locations", "1.0, 2.0");
  prm_.set ("number of cells for x, y, z directions", "1, 3");
}

TEST_F (ComputingDataTest, 2DFundamentalDataTest) {
  dealii::parallel::distributed::Triangulation<2> tria_2d(MPI_COMM_WORLD);
}
