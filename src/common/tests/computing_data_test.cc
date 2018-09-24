#include "../computing_data.h"

#include <unordered_map>
#include <vector>

#include "../problem_definition.h"
#include "../../test_helpers/bart_test_helper.h"
#include "../../material/tests/mock_material_properties.h"

// Forward declaration of test helper function
namespace btest {
std::unordered_map<int, std::vector<double>> RandomIntVectorMap();
}

class XSectionsTest : public ::testing::Test {
 protected:
  using id_vector_map = std::unordered_map<int, std::vector<double>>;
  btest::MockMaterialProperties mock_material_properties;
};

class ComputingDataTestMPI
    :
    public btest::BARTParallelEnvironment {
 public:
  void SetUp (); 
  dealii::ParameterHandler prm_;
};

void ComputingDataTestMPI::SetUp() {
  // Initilize MPI
  this->MPIInit();
  // declare all entries for parameter handler
  bparams::DeclareParameters(prm_);

  prm_.set ("reflective boundary names", "xmin");
  prm_.set ("x, y, z max values of boundary locations", "1.0, 2.0");
  prm_.set ("number of cells for x, y, z directions", "1, 3");
}

TEST_F (ComputingDataTestMPI, 2DFundamentalDataTest) {
  dealii::parallel::distributed::Triangulation<2> tria_2d(MPI_COMM_WORLD);
}

TEST_F(XSectionsTest, XSectionsConstructor) {

  id_vector_map sigma_t_map = btest::RandomIntVectorMap();
  id_vector_map sigma_t_inv_map = btest::RandomIntVectorMap();
  id_vector_map q_map = btest::RandomIntVectorMap();
  id_vector_map q_per_ster_map = btest::RandomIntVectorMap();
  id_vector_map nu_sigf_map = btest::RandomIntVectorMap();
  
  EXPECT_CALL(mock_material_properties, GetSigT())
      .WillOnce(::testing::Return(sigma_t_map));
  EXPECT_CALL(mock_material_properties, GetInvSigT())
      .WillOnce(::testing::Return(sigma_t_inv_map));
  EXPECT_CALL(mock_material_properties, GetQ())
        .WillOnce(::testing::Return(q_map));
  EXPECT_CALL(mock_material_properties, GetQPerSter())
        .WillOnce(::testing::Return(q_per_ster_map));
  EXPECT_CALL(mock_material_properties, GetNuSigF())
        .WillOnce(::testing::Return(nu_sigf_map));
  
  XSections test_xsections(mock_material_properties);

  EXPECT_EQ(test_xsections.sigt, sigma_t_map);
  EXPECT_EQ(test_xsections.inv_sigt, sigma_t_inv_map);
  EXPECT_EQ(test_xsections.q, q_map);
  EXPECT_EQ(test_xsections.q_per_ster, q_per_ster_map);
  EXPECT_EQ(test_xsections.nu_sigf, nu_sigf_map);  
}
