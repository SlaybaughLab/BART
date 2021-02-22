#include "data/cross_sections/cross_sections.hpp"

#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>
#include <deal.II/lac/full_matrix.h>

#include "problem/parameters_dealii_handler.h"
#include "material/tests/material_mock.hpp"
#include "test_helpers/bart_test_helper.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace {

using namespace bart;

class CrossSectionsTest : public ::testing::Test {
 protected:
  using id_vector_map = std::unordered_map<int, std::vector<double>>;
  using id_matrix_map = std::unordered_map<int, dealii::FullMatrix<double>>;
  ::testing::NiceMock<material::MaterialMock> mock_material_properties;
};

class CrossSectionsTestConstructor : public CrossSectionsTest {
 protected:
  void SetUp() override;
  id_vector_map diffusion_coef_map = test_helpers::RandomIntVectorMap();
  id_vector_map sigma_t_map = test_helpers::RandomIntVectorMap();
  id_vector_map sigma_t_inv_map = test_helpers::RandomIntVectorMap();
  id_vector_map q_map = test_helpers::RandomIntVectorMap();
  id_vector_map q_per_ster_map = test_helpers::RandomIntVectorMap();
  id_vector_map nu_sigf_map = test_helpers::RandomIntVectorMap();
  id_matrix_map sig_s_map = test_helpers::RandomIntMatrixMap();
  id_matrix_map sig_s_per_ster_map = test_helpers::RandomIntMatrixMap();
  id_matrix_map chi_nu_sig_f_map = test_helpers::RandomIntMatrixMap();
  id_matrix_map chi_nu_sig_f_per_ster_map = test_helpers::RandomIntMatrixMap();
  std::unordered_map<int, bool> fissile_id_map{{1, true}, {2, false}};
};

void CrossSectionsTestConstructor::SetUp() {
  ON_CALL(mock_material_properties, GetDiffusionCoef())
      .WillByDefault(::testing::Return(diffusion_coef_map));
  ON_CALL(mock_material_properties, GetSigT())
      .WillByDefault(::testing::Return(sigma_t_map));
  ON_CALL(mock_material_properties, GetInvSigT())
      .WillByDefault(::testing::Return(sigma_t_inv_map));
  ON_CALL(mock_material_properties, GetQ())
      .WillByDefault(::testing::Return(q_map));
  ON_CALL(mock_material_properties, GetQPerSter())
      .WillByDefault(::testing::Return(q_per_ster_map));
  ON_CALL(mock_material_properties, GetNuSigF())
      .WillByDefault(::testing::Return(nu_sigf_map));
  ON_CALL(mock_material_properties, GetSigS())
      .WillByDefault(::testing::Return(sig_s_map));
  ON_CALL(mock_material_properties, GetSigSPerSter())
      .WillByDefault(::testing::Return(sig_s_per_ster_map));
  ON_CALL(mock_material_properties, GetChiNuSigF())
      .WillByDefault(::testing::Return(chi_nu_sig_f_map));
  ON_CALL(mock_material_properties, GetChiNuSigFPerSter())
      .WillByDefault(::testing::Return(chi_nu_sig_f_per_ster_map));
  ON_CALL(mock_material_properties, GetFissileIDMap())
      .WillByDefault(::testing::Return(fissile_id_map));
}

TEST_F(CrossSectionsTestConstructor, CrossSectionsConstructor) {
  bart::data::cross_sections::CrossSections test_xsections(mock_material_properties);
  EXPECT_EQ(test_xsections.diffusion_coef, diffusion_coef_map);
  EXPECT_EQ(test_xsections.sigma_t, sigma_t_map);
  EXPECT_EQ(test_xsections.inverse_sigma_t, sigma_t_inv_map);
  EXPECT_EQ(test_xsections.q, q_map);
  EXPECT_EQ(test_xsections.q_per_ster, q_per_ster_map);
  EXPECT_EQ(test_xsections.nu_sigma_f, nu_sigf_map);
  EXPECT_EQ(test_xsections.sigma_s, sig_s_map);
  EXPECT_EQ(test_xsections.sigma_s_per_ster, sig_s_per_ster_map);
  EXPECT_EQ(test_xsections.fiss_transfer, chi_nu_sig_f_map);
  EXPECT_EQ(test_xsections.fiss_transfer_per_ster, chi_nu_sig_f_per_ster_map);
  EXPECT_EQ(test_xsections.is_material_fissile, fissile_id_map);

}

} // namespace