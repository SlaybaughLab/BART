#include "data/material/material_protobuf.hpp"

#include <sstream>
#include <exception>

#include <google/protobuf/text_format.h>

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace {

using namespace bart;

class MaterialProtobufTest : public ::testing::Test {
 protected:
  using MaterialProtobuf = data::material::MaterialProtobuf;
  MaterialProtobufTest();

  //! maps from material id to data for that material
  std::unordered_map<int, bool> test_is_fissile_;
  std::unordered_map<int, std::vector<double>> test_energy_groups_;
  std::unordered_map<int, std::vector<double>> test_sigma_t_;
  std::unordered_map<int, std::vector<double>> test_sigma_a_;
  std::unordered_map<int, std::vector<double>> test_nu_sig_f_;
  std::unordered_map<int, std::vector<double>> test_kappa_sig_f_;
  std::unordered_map<int, std::vector<double>> test_chi_;
  std::unordered_map<int, std::vector<double>> test_sigma_s_;

  //! working Material objects used containing above data
  std::unordered_map<int, Material> test_materials_;

  // uninitialized empty maps
  std::unordered_map<int, bool> null_bool_map_;
  std::unordered_map<int, std::vector<double>> null_vector_map_;
  std::unordered_map<int, dealii::FullMatrix<double>> null_matrix_map_;

  //! extracts the message from a dealii exception as a string and trims leading and trailing whitespace
  static std::string GetMessage(const dealii::ExceptionBase& e);

  /*!
    returns a non-const pointer to the first vector property of the material with given ID
    returns nullptr if vector property with given ID is absent
  */
  static Material_VectorProperty* GetPropertyPtr(Material& material, Material::VectorId property_id);

  /*!
    returns a non-const pointer to the first matrix property of the material with given ID
    returns nullptr if matrix property with given ID is absent
  */
  static Material_MatrixProperty* GetPropertyPtr(Material& material, Material::MatrixId property_id);

  //! compare the contents of two maps with allowance for floating point error
  static void ExpectApproxEqual(const std::unordered_map<int, dealii::FullMatrix<double>>& lhs,
                                const std::unordered_map<int, dealii::FullMatrix<double>>& rhs);

  //! compare the contents of two maps with allowance for floating point error
  static void ExpectApproxEqual(const std::unordered_map<int, std::vector<double>>& lhs,
                                const std::unordered_map<int, std::vector<double>>& rhs);
};

MaterialProtobufTest::MaterialProtobufTest() {
  /*
    test_materials_ constructed here are identical to the materials in
    BART/src/material/tests/data/serialized/ and BART/src/material/tests/data/readable/
    This can be verified by reading in the materials from those files and using
    google::protobuf::util::MessageDifferencer from <google/protobuf/util/message_differencer.h>
  */
  const std::unordered_set<int> material_ids = {1, 2, 10, 11};
  for (const int& key : material_ids) {
    test_energy_groups_[key] = {20000000.0, 1353000.0, 9119.0, 3.928, 0.6251, 0.1457, 0.05692, 0};
  }
  test_is_fissile_ = {{1, false}, {2, false}, {10, true}, {11, true}};
  test_sigma_t_[1] = {0.12417, 0.29921, 0.57997, 1.0581, 1.3203, 1.6301, 2.2847};
  test_sigma_t_[2] = {0.075384, 0.24872, 0.42163, 0.53183, 0.90849, 1.3205, 2.3163};
  test_sigma_t_[10] = {0.11111, 0.28863, 0.45098, 0.45889, 0.66863, 0.95402, 1.6043};
  test_sigma_t_[11] = {0.11113, 0.28844, 0.45382, 0.46398, 0.68795, 0.98919, 1.6809};
  test_sigma_a_[1] = {0.002142, 0.0073171, 0.18247, 0.62264, 0.68385, 0.69803, 0.6959};
  test_sigma_a_[2] = {0.00049929, 1.5705e-05, 0.00082969, 0.0054649, 0.013821, 0.021644, 0.039881};
  test_sigma_a_[10] = {0.0047082, 0.0019361, 0.023985, 0.014122, 0.043427, 0.070252, 0.1325};
  test_sigma_a_[11] = {0.0047825, 0.0020899, 0.02669, 0.018674, 0.060669, 0.09879, 0.18302};
  test_nu_sig_f_[1] = {0, 0, 0, 0, 0, 0, 0};
  test_nu_sig_f_[2] = {0, 0, 0, 0, 0, 0, 0};
  test_nu_sig_f_[10] = {0.011259, 0.00068735, 0.0077368, 0.013847, 0.060142, 0.098688, 0.18912};
  test_nu_sig_f_[11] = {0.011458, 0.001054, 0.0123, 0.022601, 0.095993, 0.15886, 0.29556};
  test_kappa_sig_f_[1] = {0, 0, 0, 0, 0, 0, 0};
  test_kappa_sig_f_[2] = {0, 0, 0, 0, 0, 0, 0};
  test_kappa_sig_f_[10] = {1.3749e-13, 9.0606e-15, 1.0307e-13, 1.8447e-13, 8.0124e-13, 1.3148e-12, 2.5195e-12};
  test_kappa_sig_f_[11] = {1.3977e-13, 1.3885e-14, 1.6387e-13, 3.011e-13, 1.2789e-12, 2.1164e-12, 3.9375e-12};
  test_chi_[10] = {0.5925247816749882, 0.4071432856463152, 0.00033193267869671705, 0, 0, 0, 0};
  test_chi_[11] = {0.5925247816749882, 0.4071432856463152, 0.00033193267869671705, 0, 0, 0, 0};
  test_sigma_s_[1] =   {
    0.12334, 0.055604999999999995, 0.00023375, 0, 0, 0, 0,
    0, 0.4042, 0.043422, 0, 0, 0, 0,
    0, 0, 0.6650699999999999, 0.041651, 0.0058453, 0.0010578, 0.00065308,
    0, 0, 0, 0.54505, 0.21211999999999998, 0.029405, 0.013149000000000001,
    0, 0, 0, 0.0043207, 0.68842, 0.22864, 0.072648,
    0, 0, 0, 0, 0.14998, 0.8552799999999999, 0.27904,
    0, 0, 0, 0, 0.06477100000000001, 0.43825, 1.4101};
  test_sigma_s_[2] = {
    0.082716, 0.081963, 0.00051642, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.47143, 0.09973, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.95552, 0.11014000000000002, 0.015750999999999998, 0.0028683000000000003, 0.0017862000000000002,
    0.0, 0.0, 0.0, 0.70714, 0.35409, 0.049957999999999995, 0.022326,
    0.0, 0.0, 0.0, 0.0023861, 0.8920299999999999, 0.42845, 0.13581,
    0.0, 0.0, 0.0, 0.0, 0.21673, 1.1939, 0.43268999999999996,
    0.0, 0.0, 0.0, 0.0, 0.098237, 0.64545, 2.0233};
  test_sigma_s_[10] = {
    0.12244000000000001, 0.067156, 0.0002876, 0, 0, 0, 0,
    0, 0.4302, 0.051664999999999996, 0, 0, 0, 0,
    0, 0, 0.7420100000000001, 0.048857, 0.0068658999999999994, 0.0012425000000000001, 0.00076711,
    0, 0, 0, 0.5464899999999999, 0.19899, 0.027304000000000002, 0.012209000000000001,
    0, 0, 0, 0.0045242, 0.66623, 0.20922, 0.064648,
    0, 0, 0, 0, 0.1404, 0.80301, 0.252,
    0, 0, 0, 0, 0.058083, 0.404, 1.299};
  test_sigma_s_[11] = {
    0.12239000000000001, 0.06713, 0.0002876, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.42991999999999997, 0.051655999999999994, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.74242, 0.048867, 0.0068715, 0.0012435, 0.00076773,
    0.0, 0.0, 0.0, 0.54727, 0.19917, 0.027343, 0.012226,
    0.0, 0.0, 0.0, 0.0047185000000000005, 0.6695, 0.20942, 0.065028,
    0.0, 0.0, 0.0, 0.0, 0.14190999999999998, 0.80934, 0.25445,
    0.0, 0.0, 0.0, 0.0, 0.05942000000000001, 0.41229, 1.3222};

  test_materials_[1].set_full_name("Control Rod");
  test_materials_[1].set_abbreviation("control_rod");
  test_materials_[2].set_full_name("Reflector");
  test_materials_[2].set_abbreviation("reflector");
  test_materials_[10].set_full_name("UO2 2.0% fuel cell");
  test_materials_[10].set_abbreviation("uo2_20");
  test_materials_[11].set_full_name("UO2 3.3% fuel cell");
  test_materials_[11].set_abbreviation("uo2_33");

  for (const int& key : material_ids) {
    test_materials_[key].set_is_fissionable(test_is_fissile_[key]);
    test_materials_[key].set_number_of_groups(7);
    test_materials_[key].set_thermal_groups(4);
    test_materials_[key].set_id(test_materials_[key].abbreviation());
    Material_VectorProperty* vec_prop_ptr;
    vec_prop_ptr = test_materials_[key].add_vector_property();
    vec_prop_ptr->set_id(Material::ENERGY_GROUPS);
    for (const double& val : test_energy_groups_[key]) {
      vec_prop_ptr->add_value(val);
    }
    vec_prop_ptr = test_materials_[key].add_vector_property();
    vec_prop_ptr->set_id(Material::SIGMA_T);
    for (const double& val : test_sigma_t_[key]) {
      vec_prop_ptr->add_value(val);
    }
    vec_prop_ptr = test_materials_[key].add_vector_property();
    vec_prop_ptr->set_id(Material::SIGMA_A);
    for (const double& val : test_sigma_a_[key]) {
      vec_prop_ptr->add_value(val);
    }
    if (test_chi_.count(key) > 0) {
      vec_prop_ptr = test_materials_[key].add_vector_property();
      vec_prop_ptr->set_id(Material::CHI);
      for (const double& val : test_chi_[key]) {
        vec_prop_ptr->add_value(val);
      }
    }
    vec_prop_ptr = test_materials_[key].add_vector_property();
    vec_prop_ptr->set_id(Material::NU_SIG_F);
    for (const double& val : test_nu_sig_f_[key]) {
      vec_prop_ptr->add_value(val);
    }
    vec_prop_ptr = test_materials_[key].add_vector_property();
    vec_prop_ptr->set_id(Material::KAPPA_SIG_F);
    for (const double& val : test_kappa_sig_f_[key]) {
      vec_prop_ptr->add_value(val);
    }
    Material_MatrixProperty* matrix_prop_ptr = test_materials_[key].add_matrix_property();
    matrix_prop_ptr->set_id(Material::SIGMA_S);
    matrix_prop_ptr->set_rows(7);
    matrix_prop_ptr->set_cols(7);
    for (const double& val : test_sigma_s_[key]) {
      matrix_prop_ptr->add_value(val);
    }
  }
}

std::string MaterialProtobufTest::GetMessage(const dealii::ExceptionBase& e) {
  std::stringstream ss;
  e.print_info(ss);
  std::string msg = ss.str();

  const std::string whitespace_chars = " \t\r\n\v\f"; // whitespace to be trimmed from start and end
  msg.erase(0, msg.find_first_not_of(whitespace_chars));
  msg.erase(msg.find_last_not_of(whitespace_chars) + 1);

  return msg;
}

Material_VectorProperty* MaterialProtobufTest::GetPropertyPtr(Material& material, Material::VectorId property_id) {
  for (int i = 0; i < material.vector_property_size(); ++i) {
    if (material.vector_property(i).id() == property_id) {
      return material.mutable_vector_property(i);
    }
  }
  return nullptr;
}

Material_MatrixProperty* MaterialProtobufTest::GetPropertyPtr(Material& material, Material::MatrixId property_id) {
  for (int i = 0; i < material.matrix_property_size(); ++i) {
    if (material.matrix_property(i).id() == property_id) {
      return material.mutable_matrix_property(i);
    }
  }
  return nullptr;
}

void MaterialProtobufTest::ExpectApproxEqual(const std::unordered_map<int, dealii::FullMatrix<double>>& lhs,
                                               const std::unordered_map<int, dealii::FullMatrix<double>>& rhs) {
  /*
    All computed values here require only a one step mathematical operation: chi*nu_sig_f, 1/sig_s, and dividing by 4pi,
    so there isn't much that can be different between implementations.
    EXPECT_DOUBLE_EQ allows error within four ULPs, but values really shouldn't be off by more than one ULP. 
  */
  std::unordered_set<int> lhs_key_set;
  std::unordered_set<int> rhs_key_set;
  for (const auto& l_pair : lhs) {
    lhs_key_set.insert(l_pair.first);
  }
  for (const auto& r_pair : rhs) {
    rhs_key_set.insert(r_pair.first);
  }

  ASSERT_EQ(lhs_key_set, rhs_key_set);
  for (const int& key : lhs_key_set) {
    SCOPED_TRACE("key = " + std::to_string(key));
    const auto& lhs_matrix = lhs.at(key);
    const auto& rhs_matrix = rhs.at(key);
    ASSERT_EQ(lhs_matrix.m(), rhs_matrix.m());
    ASSERT_EQ(lhs_matrix.n(), rhs_matrix.n());
    auto lhs_iter = lhs_matrix.begin();
    auto rhs_iter = rhs_matrix.begin();
    while (lhs_iter != lhs_matrix.end()) {
      EXPECT_DOUBLE_EQ(lhs_iter->value(), rhs_iter->value());
      ++lhs_iter;
      ++rhs_iter;
    }
  }
}

void MaterialProtobufTest::ExpectApproxEqual(const std::unordered_map<int, std::vector<double>>& lhs,
                                               const std::unordered_map<int, std::vector<double>>& rhs) {
  std::unordered_set<int> lhs_key_set;
  std::unordered_set<int> rhs_key_set;
  for (const auto& l_pair : lhs) {
    lhs_key_set.insert(l_pair.first);
  }
  for (const auto& r_pair : rhs) {
    rhs_key_set.insert(r_pair.first);
  }

  ASSERT_EQ(lhs_key_set, rhs_key_set);
  for (const int& key : lhs_key_set) {
    SCOPED_TRACE("key = " + std::to_string(key));
    const auto& lhs_vector = lhs.at(key);
    const auto& rhs_vector = rhs.at(key);
    ASSERT_EQ(lhs_vector.size(), rhs_vector.size());
    auto lhs_iter = lhs_vector.cbegin();
    auto rhs_iter = rhs_vector.cbegin();
    while (lhs_iter != lhs_vector.cend()) {
      EXPECT_DOUBLE_EQ(*lhs_iter, *rhs_iter);
      ++lhs_iter;
      ++rhs_iter;
    }
  }
}

/*
  Many of these tests test that a certain exception is thrown and that it prints the right message.
  This is done using the following format:

    EXPECT_THROW({
      try {
        // code that is supposed to throw an exception
      }
      catch (const ExpectedExceptionType& e) {
        EXPECT_EQ("expected exception message", GetMessage(e)); // assert that the message is correct
        throw; // rethrow the exception for the EXPECT_THROW
      }
    }, ExpectedExceptionType);

  If an exception with the wrong message is thrown, the EXPECT_EQ will fail,
  and if no exception of the correct type is thrown, the EXPECT_THROW will fail.
*/

// tests for static MaterialProtobuf::CheckValid, named after exceptions

TEST_F(MaterialProtobufTest, MultipleDefinition) {
  Material uo2_20;
  std::string uo2_20_msg = "Material full_name|abbreviation|id = \"UO2 2.0% fuel cell\"|\"uo2_20\"|\"uo2_20\"";

  EXPECT_THROW({
    try {
      uo2_20.CopyFrom(test_materials_.at(10));
      Material_VectorProperty* extra_vec_prop = uo2_20.add_vector_property();
      extra_vec_prop->set_id(Material::ENERGY_GROUPS);
      const std::vector<double> energy_groups({1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 10, 1, 0});
      for (const double& val : energy_groups) {
        extra_vec_prop->add_value(val);
      }
      MaterialProtobuf::CheckValid(uo2_20);
    }
    catch (const MaterialProtobuf::MultipleDefinition& e) {
      std::string expected = uo2_20_msg + " has multiple definitions of ENERGY_GROUPS, which is not allowed.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::MultipleDefinition);

  EXPECT_THROW({
    try {
      uo2_20.CopyFrom(test_materials_.at(10));
      Material_MatrixProperty* extra_mat_prop = uo2_20.add_matrix_property();
      extra_mat_prop->set_id(Material::SIGMA_S);
      for (unsigned int i = 0; i < uo2_20.number_of_groups()*uo2_20.number_of_groups(); ++i) {
        extra_mat_prop->add_value(0);
      }
      MaterialProtobuf::CheckValid(uo2_20, false);
    }
    catch (const MaterialProtobuf::MultipleDefinition& e) {
      std::string expected = uo2_20_msg + " has multiple definitions of SIGMA_S, which is not allowed.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::MultipleDefinition);

  EXPECT_THROW({
    try {
      uo2_20.CopyFrom(test_materials_.at(10));
      Material_VectorProperty* extra_vec_prop = uo2_20.add_vector_property();
      extra_vec_prop->set_id(Material::CHI);
      const std::vector<double> chi({0.8, 0.2, 0, 0, 0, 0, 0});
      for (const double& val : chi) {
        extra_vec_prop->add_value(val);
      }
      MaterialProtobuf::CheckValid(uo2_20, true);
    }
    catch (const MaterialProtobuf::MultipleDefinition& e) {
      std::string expected = uo2_20_msg + " has multiple definitions of CHI, which is not allowed.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::MultipleDefinition);
}

TEST_F(MaterialProtobufTest, MissingProperty) {
  Material uo2_33;
  std::string uo2_33_msg = "Material full_name|abbreviation|id = \"UO2 3.3% fuel cell\"|\"uo2_33\"|\"uo2_33\"";

  EXPECT_THROW({
    try {
      uo2_33.CopyFrom(test_materials_.at(11));
      GetPropertyPtr(uo2_33, Material::SIGMA_T)->set_id(Material::UNKNOWN_VECTOR);
      MaterialProtobuf::CheckValid(uo2_33);
    }
    catch (const MaterialProtobuf::MissingProperty& e) {
      std::string expected = uo2_33_msg + " is missing required property SIGMA_T.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::MissingProperty);

  EXPECT_THROW({
    try {
      uo2_33.CopyFrom(test_materials_.at(11));
      uo2_33.clear_matrix_property();
      MaterialProtobuf::CheckValid(uo2_33, false);
    }
    catch (const MaterialProtobuf::MissingProperty& e) {
      std::string expected = uo2_33_msg + " is missing required property SIGMA_S.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::MissingProperty);

  EXPECT_THROW({
    try {
      uo2_33.CopyFrom(test_materials_.at(11));
      GetPropertyPtr(uo2_33, Material::NU_SIG_F)->set_id(Material::Q);
      MaterialProtobuf::CheckValid(uo2_33, true);
    }
    catch (const MaterialProtobuf::MissingProperty& e) {
      std::string expected = uo2_33_msg + " is missing required property NU_SIG_F.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::MissingProperty);
}

TEST_F(MaterialProtobufTest, WrongNumberOfValues) {
  Material control_rod;
  Material uo2_33;
  std::string control_rod_msg = "In material full_name|abbreviation|id = \"Control Rod\"|\"control_rod\"|\"control_rod\"";

  EXPECT_THROW({
    try {
      uo2_33.CopyFrom(test_materials_.at(11));
      GetPropertyPtr(uo2_33, Material::NU_SIG_F)->clear_value();
      MaterialProtobuf::CheckValid(uo2_33, true);
    }
    catch (const MaterialProtobuf::WrongNumberOfValues& e) {
      std::string expected = "In material full_name|abbreviation|id = \"UO2 3.3% fuel cell\"|\"uo2_33\"|\"uo2_33\"";
      expected += ", the number of values under the label NU_SIG_F (0) is inconsistent with number_of_groups = 7.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongNumberOfValues);

  EXPECT_THROW({
    try {
      uo2_33.CopyFrom(test_materials_.at(11));
      uo2_33.set_number_of_groups(8);
      MaterialProtobuf::CheckValid(uo2_33);
    }
    catch (const MaterialProtobuf::WrongNumberOfValues& e) {
      EXPECT_EQ('8', *(GetMessage(e).rbegin()+1));
      throw;
    }
  }, MaterialProtobuf::WrongNumberOfValues);

  EXPECT_THROW({
    try {
      control_rod.CopyFrom(test_materials_.at(1));
      GetPropertyPtr(control_rod, Material::SIGMA_S)->add_value(0.1);
      MaterialProtobuf::CheckValid(control_rod, false);
    }
    catch (const MaterialProtobuf::WrongNumberOfValues& e) {
      std::string expected = control_rod_msg;
      expected += ", the number of values under the label SIGMA_S (50) is inconsistent with number_of_groups = 7.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongNumberOfValues);

  EXPECT_THROW({
    try {
      control_rod.CopyFrom(test_materials_.at(1));
      GetPropertyPtr(control_rod, Material::SIGMA_S)->clear_value();
      for (int i = 0; i < 36; ++i) {
        GetPropertyPtr(control_rod, Material::SIGMA_S)->add_value(0.1);
      }
      MaterialProtobuf::CheckValid(control_rod, false);
    }
    catch (const MaterialProtobuf::WrongNumberOfValues& e) {
      std::string expected = control_rod_msg;
      expected += ", the number of values under the label SIGMA_S (36) is inconsistent with number_of_groups = 7.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongNumberOfValues);

  EXPECT_THROW({
    try {
      control_rod.CopyFrom(test_materials_.at(1));
      GetPropertyPtr(control_rod, Material::ENERGY_GROUPS)->clear_value();
      const std::vector<double> energy_groups({1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 0});
      for (const double& val : energy_groups) {
        GetPropertyPtr(control_rod, Material::ENERGY_GROUPS)->add_value(val);
      }
      MaterialProtobuf::CheckValid(control_rod, false);
    }
    catch (const MaterialProtobuf::WrongNumberOfValues& e) {
      std::string expected = control_rod_msg;
      expected += ", the number of values under the label ENERGY_GROUPS (7) is inconsistent with number_of_groups = 7.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongNumberOfValues);

  // if the material contains Q data, it must be checked
  EXPECT_THROW({
    try {
      control_rod.CopyFrom(test_materials_.at(1));
      Material_VectorProperty* q_ptr = control_rod.add_vector_property();
      q_ptr->set_id(Material::Q);
      const std::vector<double> fixed_source({1, 2, 3, 4, 5, 6, 7, 8});
      for (const double& val : fixed_source) {
        q_ptr->add_value(val);
      }
      MaterialProtobuf::CheckValid(control_rod, false);
    }
    catch (const MaterialProtobuf::WrongNumberOfValues& e) {
      std::string expected = control_rod_msg;
      expected += ", the number of values under the label Q (8) is inconsistent with number_of_groups = 7.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongNumberOfValues);
}

TEST_F(MaterialProtobufTest, WrongSign) {
  std::string control_rod_msg = "Material full_name|abbreviation|id = \"Control Rod\"|\"control_rod\"|\"control_rod\"";
  Material control_rod;
  control_rod.CopyFrom(test_materials_.at(1));
  GetPropertyPtr(control_rod, Material::SIGMA_T)->set_value(6, -1e-16);

  EXPECT_THROW({
    try {
      MaterialProtobuf::CheckValid(control_rod);
    }
    catch (const MaterialProtobuf::WrongSign& e) {
      std::string expected = control_rod_msg + " contains a value of the wrong sign in SIGMA_T. (-1e-16 at index 6)";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongSign);

  GetPropertyPtr(control_rod, Material::SIGMA_T)->set_value(6, 2.2847);
  GetPropertyPtr(control_rod, Material::SIGMA_S)->set_value(30, -85223);

  EXPECT_THROW({
    try {
      MaterialProtobuf::CheckValid(control_rod, false);
    }
    catch (const MaterialProtobuf::WrongSign& e) {
      std::string expected = control_rod_msg + " contains a value of the wrong sign in SIGMA_S. (-85223 at index 30)";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongSign);

  Material uo2_33;
  uo2_33.CopyFrom(test_materials_.at(11));
  GetPropertyPtr(uo2_33, Material::CHI)->set_value(4, -1e-16);

  EXPECT_THROW({
    try {
      MaterialProtobuf::CheckValid(uo2_33, true);
    }
    catch (const MaterialProtobuf::WrongSign& e) {
      std::string expected = "Material full_name|abbreviation|id = \"UO2 3.3% fuel cell\"|\"uo2_33\"|\"uo2_33\"";
      expected += " contains a value of the wrong sign in CHI. (-1e-16 at index 4)";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongSign);

  EXPECT_THROW({
    try {
      control_rod.CopyFrom(test_materials_.at(1));
      Material_VectorProperty* q_ptr = control_rod.add_vector_property();
      q_ptr->set_id(Material::Q);
      const std::vector<double> fixed_source({1, 2, 3, 4, 5, -6, 7});
      for (const double& val : fixed_source) {
        q_ptr->add_value(val);
      }
      MaterialProtobuf::CheckValid(control_rod, false);
    }
    catch (const MaterialProtobuf::WrongSign& e) {
      std::string expected = control_rod_msg + " contains a value of the wrong sign in Q. (-6 at index 5)";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongSign);
}

TEST_F(MaterialProtobufTest, EnergyGroupBoundariesNotSorted) {
  Material control_rod;
  control_rod.CopyFrom(test_materials_.at(1));
  std::string control_rod_msg = "Energy group boundary vector for material full_name|abbreviation|id =";
  control_rod_msg += " \"Control Rod\"|\"control_rod\"|\"control_rod\"";

  EXPECT_THROW({
    try {
      GetPropertyPtr(control_rod, Material::ENERGY_GROUPS)->clear_value();
      const std::vector<double> energy_groups({1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 0, 10});
      for (const double& val : energy_groups) {
        GetPropertyPtr(control_rod, Material::ENERGY_GROUPS)->add_value(val);
      }
      MaterialProtobuf::CheckValid(control_rod);
    }
    catch (const MaterialProtobuf::EnergyGroupBoundariesNotSorted& e) {
      std::string expected = control_rod_msg + " contains {100, 0, 10}, which are out of order.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::EnergyGroupBoundariesNotSorted);

  EXPECT_THROW({
    try {
      GetPropertyPtr(control_rod, Material::ENERGY_GROUPS)->clear_value();
      const std::vector<double> energy_groups({0, 10, 1e2, 1e4, 1e3, 1e5, 1e6, 1e7});
      for (const double& val : energy_groups) {
        GetPropertyPtr(control_rod, Material::ENERGY_GROUPS)->add_value(val);
      }
      MaterialProtobuf::CheckValid(control_rod);
    }
    catch (const MaterialProtobuf::EnergyGroupBoundariesNotSorted& e) {
      std::string expected = control_rod_msg + " contains {100, 10000, 1000}, which are out of order.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::EnergyGroupBoundariesNotSorted);
}

TEST_F(MaterialProtobufTest, ChiDoesNotSumToOne) {
  Material uo2_33;
  uo2_33.CopyFrom(test_materials_.at(11));

  EXPECT_THROW({
    try {
      GetPropertyPtr(uo2_33, Material::CHI)->clear_value();
      const std::vector<double> chi({0.5, 0.5, 0, 0, 0.5, 0.25, 0.25});
      for (const double& val : chi) {
        GetPropertyPtr(uo2_33, Material::CHI)->add_value(val);
      }
      MaterialProtobuf::CheckValid(uo2_33, true);
    }
    catch (const MaterialProtobuf::ChiDoesNotSumToOne& e) {
      std::string expected = "In material full_name|abbreviation|id = \"UO2 3.3% fuel cell\"|\"uo2_33\"|\"uo2_33\"";
      expected += ", the sum of values in the Chi vector is 2. It must be normalized to 1.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::ChiDoesNotSumToOne);

  GetPropertyPtr(uo2_33, Material::CHI)->clear_value();
  const std::vector<double> chi({0.5925247816749882 - 1e-15, 0.4071432856463152, 0.00033193267869671705, 0, 0, 0, 0});
  for (const double& val : chi) {
    GetPropertyPtr(uo2_33, Material::CHI)->add_value(val);
  }

  EXPECT_THROW(MaterialProtobuf::CheckValid(uo2_33, true), MaterialProtobuf::ChiDoesNotSumToOne);
}

TEST_F(MaterialProtobufTest, CheckValidNoThrow) {
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(10), false));
  GetPropertyPtr(test_materials_.at(10), Material::KAPPA_SIG_F)->clear_value();
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(10), true));
  GetPropertyPtr(test_materials_.at(10), Material::SIGMA_A)->set_value(1, -1);
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(10), true));
  GetPropertyPtr(test_materials_.at(10), Material::SIGMA_A)->clear_value();
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(10), true));
  GetPropertyPtr(test_materials_.at(10), Material::KAPPA_SIG_F)->set_id(Material::UNKNOWN_VECTOR);
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(10), true));
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(1), false));
  
  Material_VectorProperty* unknown_vector_ptr = test_materials_.at(1).add_vector_property();
  unknown_vector_ptr->set_id(Material::UNKNOWN_VECTOR);
  unknown_vector_ptr = test_materials_.at(1).add_vector_property();
  unknown_vector_ptr->set_id(Material::UNKNOWN_VECTOR);
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(1), false));
  GetPropertyPtr(test_materials_.at(1), Material::KAPPA_SIG_F)->set_value(4, -510);
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(1), false));
  GetPropertyPtr(test_materials_.at(11), Material::KAPPA_SIG_F)->add_value(-44);
  EXPECT_NO_THROW(MaterialProtobuf::CheckValid(test_materials_.at(11), true));
}

TEST_F(MaterialProtobufTest, StaticCheckConsistent) {
  Material uo2_20;
  uo2_20.CopyFrom(test_materials_.at(10));
  Material uo2_33;
  uo2_33.CopyFrom(test_materials_.at(11));
  std::unordered_map<int, Material> map = {{3, uo2_20}, {12, uo2_33}};
  std::string mat_3 = "material number|full_name|abbreviation|id = 3|\"UO2 2.0% fuel cell\"|\"uo2_20\"|\"uo2_20\"";
  std::string mat_12 = "material number|full_name|abbreviation|id = 12|\"UO2 3.3% fuel cell\"|\"uo2_33\"|\"uo2_33\"";

  map.at(3).set_number_of_groups(6);

  EXPECT_THROW({
    try {
      MaterialProtobuf::CheckConsistent(map);
    }
    catch (const MaterialProtobuf::NumberOfGroupsMismatch& e) {
      std::string expected_1 = "The number_of_groups in " + mat_3 + " and " + mat_12 + " are different.";
      std::string expected_2 = "The number_of_groups in " + mat_12 + " and " + mat_3 + " are different.";
      EXPECT_TRUE(expected_1 == GetMessage(e) || expected_2 == GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::NumberOfGroupsMismatch);

  map.at(3).set_number_of_groups(7);
  GetPropertyPtr(map.at(12), Material::ENERGY_GROUPS)->set_value(2, 9118);

  EXPECT_THROW({
    try {
      MaterialProtobuf::CheckConsistent(map);
    }
    catch (const MaterialProtobuf::NumberOfGroupsMismatch& e) {
      std::string expected_1 = "The energy groups boundaries in " + mat_3 + " and " + mat_12 + " differ at index 2.";
      std::string expected_2 = "The energy groups boundaries in " + mat_12 + " and " + mat_3 + " differ at index 2.";
      EXPECT_TRUE(expected_1 == GetMessage(e) || expected_2 == GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::EnergyGroupsMismatch);

  GetPropertyPtr(map.at(12), Material::ENERGY_GROUPS)->set_value(2, 9119);
  GetPropertyPtr(map.at(12), Material::ENERGY_GROUPS)->add_value(0);
  
  EXPECT_THROW({
    try {
      MaterialProtobuf::CheckConsistent(map);
    }
    catch (const MaterialProtobuf::NumberOfGroupsMismatch& e) {
      std::string expected_1 = "The energy groups boundaries in " + mat_3 + " and " + mat_12 + " differ at index 9.";
      std::string expected_2 = "The energy groups boundaries in " + mat_12 + " and " + mat_3 + " differ at index 9.";
      EXPECT_TRUE(expected_1 == GetMessage(e) || expected_2 == GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::EnergyGroupsMismatch);
}

TEST_F(MaterialProtobufTest, StaticCheckConsistentNoThrow) {
  EXPECT_NO_THROW(MaterialProtobuf::CheckConsistent(test_materials_));
  EXPECT_NO_THROW(MaterialProtobuf::CheckConsistent({}));
  EXPECT_NO_THROW(MaterialProtobuf::CheckConsistent({{1, test_materials_.at(1)}}));
  test_materials_.at(1).clear_matrix_property();
  EXPECT_NO_THROW(MaterialProtobuf::CheckConsistent(test_materials_));
  EXPECT_NO_THROW(MaterialProtobuf::CheckConsistent({{1, Material()}}));
  test_materials_.at(1).set_number_of_groups(5);
  test_materials_.at(2).set_number_of_groups(5);
  EXPECT_NO_THROW(MaterialProtobuf::CheckConsistent({{1, test_materials_.at(1)}, {2, test_materials_.at(2)}}));
  GetPropertyPtr(test_materials_.at(1), Material::ENERGY_GROUPS)->clear_value();
  GetPropertyPtr(test_materials_.at(2), Material::ENERGY_GROUPS)->clear_value();
  EXPECT_NO_THROW(MaterialProtobuf::CheckConsistent({{1, test_materials_.at(1)}, {2, test_materials_.at(2)}}));
  test_materials_.at(1).clear_vector_property();
  test_materials_.at(2).clear_vector_property();
  EXPECT_NO_THROW(MaterialProtobuf::CheckConsistent({{1, test_materials_.at(1)}, {2, test_materials_.at(2)}}));
}

// tests for exceptions thrown from MaterialProtobuf Material map constructor

TEST_F(MaterialProtobufTest, NoFissileIDs) {
  test_materials_.erase(10);
  test_materials_.erase(11);
  EXPECT_THROW(MaterialProtobuf mp(test_materials_, true, false, 7, 4), MaterialProtobuf::NoFissileIDs);
}

TEST_F(MaterialProtobufTest, WrongNumberOfMaterials) {
  EXPECT_THROW({
    try {
      MaterialProtobuf mp(test_materials_, true, false, 7, 3);
    }
    catch (const MaterialProtobuf::WrongNumberOfMaterials& e) {
      std::string expected = "The actual number of materials read in (4) does not match the number of materials specified (3).";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongNumberOfMaterials);
}

TEST_F(MaterialProtobufTest, WrongNumberOfGroups) {
  Material uo2_20;
  uo2_20.CopyFrom(test_materials_.at(10));
  std::unordered_map<int, Material> map = {{10, uo2_20}};
  EXPECT_THROW({
    try {
      MaterialProtobuf mp(map, false, false, 8, 1);
    }
    catch (const MaterialProtobuf::WrongNumberOfGroups& e) {
      std::string expected = "The number_of_groups in material";
      expected += " number|full_name|abbreviation|id = 10|\"UO2 2.0% fuel cell\"|\"uo2_20\"|\"uo2_20\"";
      expected += " does not match the number of groups in MaterialProtobuf. (7 != 8)";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongNumberOfGroups);
}

TEST_F(MaterialProtobufTest, InvalidMaterialInConstructor) {
  std::string uo2_33_msg = "Material number|full_name|abbreviation|id = 11|\"UO2 3.3% fuel cell\"|\"uo2_33\"|\"uo2_33\"";
  std::string reflector_msg = "Material number|full_name|abbreviation|id = 2|\"Reflector\"|\"reflector\"|\"reflector\"";
  GetPropertyPtr(test_materials_.at(11), Material::NU_SIG_F)->set_id(Material::UNKNOWN_VECTOR);

  EXPECT_THROW({
    try {
      MaterialProtobuf mp(test_materials_, true, false, 7, 4);
    }
    catch (const MaterialProtobuf::MissingProperty& e) {
      std::string expected = uo2_33_msg + " is missing required property NU_SIG_F.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::MissingProperty);

  GetPropertyPtr(test_materials_.at(11), Material::UNKNOWN_VECTOR)->set_id(Material::NU_SIG_F);
  GetPropertyPtr(test_materials_.at(2), Material::SIGMA_S)->set_value(48, -1);

  EXPECT_THROW({
    try {
      MaterialProtobuf mp(test_materials_, true, false, 7, 4);
    }
    catch (const MaterialProtobuf::WrongSign& e) {
      std::string expected = reflector_msg + " contains a value of the wrong sign in SIGMA_S. (-1 at index 48)";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::WrongSign);
}

TEST_F(MaterialProtobufTest, InconsistentMaterialsInConstructor) {
  Material reflector;
  reflector.CopyFrom(test_materials_.at(2));
  Material uo2_33;
  uo2_33.CopyFrom(test_materials_.at(11));
  std::unordered_map<int, Material> map = {{0, reflector}, {1, uo2_33}};
  std::string mat_0 = "material number|full_name|abbreviation|id = 0|\"Reflector\"|\"reflector\"|\"reflector\"";
  std::string mat_1 = "material number|full_name|abbreviation|id = 1|\"UO2 3.3% fuel cell\"|\"uo2_33\"|\"uo2_33\"";
  GetPropertyPtr(map.at(0), Material::ENERGY_GROUPS)->set_value(0, 3e7);

  EXPECT_THROW({
    try {
      MaterialProtobuf mp(map, true, false, 7, 2);
    }
    catch (const MaterialProtobuf::NumberOfGroupsMismatch& e) {
      std::string expected_1 = "The energy groups boundaries in " + mat_0 + " and " + mat_1 + " differ at index 0.";
      std::string expected_2 = "The energy groups boundaries in " + mat_1 + " and " + mat_0 + " differ at index 0.";
      EXPECT_TRUE(expected_1 == GetMessage(e) || expected_2 == GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::EnergyGroupsMismatch);
}

// tests for MaterialProtobuf working correctly and returning the correct data given a map of Materials

TEST_F(MaterialProtobufTest, NullMaterialProtobuf) {
  MaterialProtobuf mp_0(std::unordered_map<int, Material>(), false, false, 0, 0);

  EXPECT_EQ(mp_0.GetFissileIDMap(), null_bool_map_);
  EXPECT_EQ(mp_0.GetSigT(), null_vector_map_);
  EXPECT_EQ(mp_0.GetInvSigT(), null_vector_map_);
  EXPECT_EQ(mp_0.GetQ(), null_vector_map_);
  EXPECT_EQ(mp_0.GetQPerSter(), null_vector_map_);
  EXPECT_EQ(mp_0.GetNuSigF(), null_vector_map_);
  EXPECT_EQ(mp_0.GetSigS(), null_matrix_map_);
  EXPECT_EQ(mp_0.GetSigSPerSter(), null_matrix_map_);
  EXPECT_EQ(mp_0.GetChiNuSigF(), null_matrix_map_);
  EXPECT_EQ(mp_0.GetChiNuSigFPerSter(), null_matrix_map_);
}

TEST_F(MaterialProtobufTest, ConstructorFromMap) {
  std::unordered_map<int, bool> correct_fissile_id_map;
  std::unordered_map<int, std::vector<double>> correct_sig_t_map;
  std::unordered_map<int, std::vector<double>> correct_inv_sig_t_map;
  std::unordered_map<int, std::vector<double>> correct_q_map;
  std::unordered_map<int, std::vector<double>> correct_q_per_ster_map;
  std::unordered_map<int, std::vector<double>> correct_nu_sig_f_map;
  std::unordered_map<int, std::vector<double>> sigma_s_values;
  std::unordered_map<int, std::vector<double>> sigma_s_per_ster_values;
  std::unordered_map<int, std::vector<double>> chi_nu_sig_f_values;
  std::unordered_map<int, std::vector<double>> chi_nu_sig_f_per_ster_values;
  std::unordered_map<int, dealii::FullMatrix<double>> correct_sigma_s_map;
  std::unordered_map<int, dealii::FullMatrix<double>> correct_sigma_s_per_ster_map;
  std::unordered_map<int, dealii::FullMatrix<double>> correct_chi_nu_sig_f_map;
  std::unordered_map<int, dealii::FullMatrix<double>> correct_chi_nu_sig_f_per_ster_map;
  
  // test with one Material, no Q data
  Material control_rod;
  control_rod.CopyFrom(test_materials_.at(1));
  std::unordered_map<int, Material> map = {{1, control_rod}};
  MaterialProtobuf mp_1(map, false, false, 7, 1);

  correct_fissile_id_map[1] = false;
  correct_sig_t_map[1] = test_sigma_t_.at(1);
  correct_inv_sig_t_map[1] =
   {8.053475074494644, 3.34213428695565, 1.7242271151956137,
    0.9450902561194594, 0.7574036203893054, 0.6134592969756456, 0.43769422681314835};
  correct_q_map[1] = {0, 0, 0, 0, 0, 0, 0};
  correct_q_per_ster_map[1] = {0, 0, 0, 0, 0, 0, 0};
  sigma_s_values[1] = test_sigma_s_.at(1);
  sigma_s_per_ster_values[1] = {
    0.009815085340477186, 0.00442490530531242, 1.8601233973865266e-05, 0, 0, 0, 0,
    0, 0.032165213998872053, 0.00345541296946814, 0, 0, 0, 0,
    0, 0, 0.052924589001063414, 0.0033144812673602665, 0.00046515419442752796, 8.417704940130345e-05, 5.197045511722751e-05,
    0, 0, 0, 0.04337370086611878, 0.016879973264326418, 0.0023399755508085912, 0.001046364173357666,
    0, 0, 0, 0.0003438303813085761, 0.0547827229616613, 0.018194593094265476, 0.005781144152870007,
    0, 0, 0, 0, 0.011935029182461232, 0.06806101986381811, 0.02220529766018124,
    0, 0, 0, 0, 0.005154312409502577, 0.03487482690501156, 0.11221219262694081};

  correct_sigma_s_map[1] = dealii::FullMatrix<double>(7, 7, sigma_s_values[1].data());
  correct_sigma_s_per_ster_map[1] = dealii::FullMatrix<double>(7, 7, sigma_s_per_ster_values[1].data());

  // the SCOPED_TRACE macro makes it so that if a test fails inside ExpectApproxEqual,
  // the line number will be included in the message
  EXPECT_EQ(mp_1.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_1.GetSigT(), correct_sig_t_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_1.GetQ(), correct_q_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1.GetQPerSter(), correct_q_per_ster_map);};
  EXPECT_EQ(mp_1.GetNuSigF(), null_vector_map_);
  EXPECT_EQ(mp_1.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  EXPECT_EQ(mp_1.GetChiNuSigF(), null_matrix_map_);
  EXPECT_EQ(mp_1.GetChiNuSigFPerSter(), null_matrix_map_);
  EXPECT_EQ(mp_1.GetSigS().at(1)(0, 1), 0.055604999999999995);

  // one material with Q data
  Material_VectorProperty* q_ptr = map[1].add_vector_property();
  q_ptr->set_id(Material::Q);
  const std::vector<double> fixed_source({0.05, 0.3, 0.5, 0.1, 0, 0, 0});
  for (const double& val : fixed_source) {
    q_ptr->add_value(val);
  }
  correct_q_map[1] = {0.05, 0.3, 0.5, 0.1, 0, 0, 0};
  correct_q_per_ster_map[1] =
   {0.0039788735772973835, 0.0238732414637843, 0.039788735772973836, 0.007957747154594767, 0, 0, 0};
  MaterialProtobuf mp_1_q(map, false, false, 7, 1);

  EXPECT_EQ(mp_1_q.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_1_q.GetSigT(), correct_sig_t_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1_q.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_1_q.GetQ(), correct_q_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1_q.GetQPerSter(), correct_q_per_ster_map);}
  EXPECT_EQ(mp_1_q.GetNuSigF(), null_vector_map_);
  EXPECT_EQ(mp_1_q.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1_q.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  EXPECT_EQ(mp_1_q.GetChiNuSigF(), null_matrix_map_);
  EXPECT_EQ(mp_1_q.GetChiNuSigFPerSter(), null_matrix_map_);
  EXPECT_DOUBLE_EQ(mp_1_q.GetQPerSter().at(1)[0], 0.0039788735772973835);

  // with two materials, one with Q data, one without
  map[10] = Material();
  map[10].CopyFrom(test_materials_.at(10));
  map[10].set_is_fissionable(false);

  correct_fissile_id_map[10] = false;
  correct_sig_t_map[10] = test_sigma_t_.at(10);
  correct_inv_sig_t_map[10] =
   {9.000090000900009, 3.464643314970724, 2.2173932325158545,
    2.1791714790036827, 1.4955954713369128, 1.0481960545900506, 0.6233248145608676};
  correct_q_map[10] = {0, 0, 0, 0, 0, 0, 0};
  correct_q_per_ster_map[10] = {0, 0, 0, 0, 0, 0, 0};
  sigma_s_values[10] = test_sigma_s_.at(10);
  sigma_s_per_ster_values[10] = {
    0.009743465616085833, 0.00534410467913966173, 2.288648081661455e-05, 0, 0, 0, 0,
    0, 0.03423422825906669, 0.004111370067421386, 0, 0, 0, 0,
    0, 0, 0.05904727966180864, 0.0038879165273203653, 0.0005463709618873221, 9.887500839583999e-05, 6.104467419761192e-05,
    0, 0, 0, 0.04348829242514494, 0.01583512106292813, 0.0021727832830905555, 0.00097156135010447526,
    0, 0, 0, 0.00036002439676817644, 0.05301689886805672, 0.016649198596843173, 0.005144524380502425,
    0, 0, 0, 0, 0.011172677005051052, 0.06390150542611144, 0.020053522829578813,
    0, 0, 0, 0, 0.0046220982798032785, 0.032149298504562863, 0.10337113553818602};

  correct_sigma_s_map[10] =
  {dealii::FullMatrix<double>(7, 7, sigma_s_values[10].data())};
  correct_sigma_s_per_ster_map[10] =
  {dealii::FullMatrix<double>(7, 7, sigma_s_per_ster_values[10].data())};

  MaterialProtobuf mp_2(map, false, false, 7, 2);

  EXPECT_EQ(mp_2.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_2.GetSigT(), correct_sig_t_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_2.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_2.GetQ(), correct_q_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_2.GetQPerSter(), correct_q_per_ster_map);}
  EXPECT_EQ(mp_2.GetNuSigF(), null_vector_map_);
  EXPECT_EQ(mp_2.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_2.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  EXPECT_EQ(mp_2.GetChiNuSigF(), null_matrix_map_);
  EXPECT_EQ(mp_2.GetChiNuSigFPerSter(), null_matrix_map_);

  // two material eigen problem

  correct_fissile_id_map[10] = true;
  correct_nu_sig_f_map[10] = test_nu_sig_f_.at(10);
  chi_nu_sig_f_values[10] = {
   0.006671236516878692, 0.004584026253091862, 3.737230029446337e-06, 0, 0, 0, 0,
   0.00040727190868430316, 0.00027984993738899475, 2.2815392670218845e-07, 0, 0, 0, 0,
   0.004584245730863049, 0.0031499861723884113, 2.5680967485407607e-06, 0, 0, 0, 0,
   0.008204690651853561, 0.005637713076344526, 4.596271801913441e-06, 0, 0, 0, 0,
   0.03563562541949714, 0.024486411485340687, 1.9963095162177958e-05, 0, 0, 0, 0,
   0.05847508565394124, 0.040180156573863555, 3.275777219522161e-05, 0, 0, 0, 0,
   0.11205828671037378, 0.07699893818143114, 6.277510819512313e-05, 0, 0, 0, 0};

  chi_nu_sig_f_per_ster_values[10] = {
   0.0005308801340982011, 0.0003647852187212948, 2.973993163289271e-07, 0.0, 0.0, 0.0, 0.0,
   3.2409668724788937e-05, 2.2269750429707967e-05, 1.815591261023963e-08, 0.0, 0.0, 0.0, 0.0,
   0.00036480268420738634, 0.0002506679350033674, 2.043626459342431e-07, 0.0, 0.0, 0.0, 0.0,
   0.0006529085368911796, 0.00044863495191702367, 3.657596885342085e-07, 0.0, 0.0, 0.0, 0.0,
   0.0028357929678420836, 0.0019485667132370648, 1.588612637237262e-06, 0.0, 0.0, 0.0, 0.0,
   0.004653299464773362, 0.0031974352664683493, 2.6067806847738836e-06, 0.0, 0.0, 0.0, 0.0,
   0.008917315122182416, 0.00612738081220102, 4.995484386191197e-06, 0.0, 0.0, 0.0, 0.0};

  correct_chi_nu_sig_f_map[10] = dealii::FullMatrix<double>(7, 7, chi_nu_sig_f_values[10].data());
  correct_chi_nu_sig_f_per_ster_map[10] = dealii::FullMatrix<double>(7, 7, chi_nu_sig_f_per_ster_values[10].data());
  
  map[10].set_is_fissionable(true);
  MaterialProtobuf mp_2_eigen(map, true, false, 7, 2);

  EXPECT_EQ(mp_2_eigen.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_2_eigen.GetSigT(), correct_sig_t_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_2_eigen.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_2_eigen.GetQ(), null_vector_map_);
  EXPECT_EQ(mp_2_eigen.GetQPerSter(), null_vector_map_);
  EXPECT_EQ(mp_2_eigen.GetNuSigF(), correct_nu_sig_f_map);
  EXPECT_EQ(mp_2_eigen.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_2_eigen.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_2_eigen.GetChiNuSigF(), correct_chi_nu_sig_f_map);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_2_eigen.GetChiNuSigFPerSter(), correct_chi_nu_sig_f_per_ster_map);}
  EXPECT_EQ(mp_2_eigen.GetNuSigF().at(10)[0], 0.011259);
  EXPECT_DOUBLE_EQ(mp_2_eigen.GetChiNuSigF().at(10)(6, 1), 0.07699893818143114);
  EXPECT_DOUBLE_EQ(mp_2_eigen.GetChiNuSigFPerSter().at(10)(5, 2), 2.6067806847738836e-06);

  // one material eigen problem
  auto correct_fissile_id_map_1_eigen = correct_fissile_id_map;
  correct_fissile_id_map_1_eigen.erase(1);
  auto correct_sig_t_map_1_eigen = correct_sig_t_map;
  correct_sig_t_map_1_eigen.erase(1);
  auto correct_inv_sig_t_map_1_eigen = correct_inv_sig_t_map;
  correct_inv_sig_t_map_1_eigen.erase(1);
  auto correct_nu_sig_f_map_1_eigen = correct_nu_sig_f_map;
  correct_nu_sig_f_map_1_eigen.erase(1);
  auto correct_sigma_s_map_1_eigen = correct_sigma_s_map;
  correct_sigma_s_map_1_eigen.erase(1);
  auto correct_sigma_s_per_ster_map_1_eigen = correct_sigma_s_per_ster_map;
  correct_sigma_s_per_ster_map_1_eigen.erase(1);
  auto correct_chi_nu_sig_f_map_1_eigen = correct_chi_nu_sig_f_map;
  correct_chi_nu_sig_f_map_1_eigen.erase(1);
  auto correct_chi_nu_sig_f_per_ster_map_1_eigen = correct_chi_nu_sig_f_per_ster_map;
  correct_chi_nu_sig_f_per_ster_map_1_eigen.erase(1);
  auto map_1_eigen = map;
  map_1_eigen.erase(1);
  MaterialProtobuf mp_1_eigen(map_1_eigen, true, false, 7, 1);

  EXPECT_EQ(mp_1_eigen.GetFissileIDMap(), correct_fissile_id_map_1_eigen);
  EXPECT_EQ(mp_1_eigen.GetSigT(), correct_sig_t_map_1_eigen);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1_eigen.GetInvSigT(), correct_inv_sig_t_map_1_eigen);}
  EXPECT_EQ(mp_1_eigen.GetQ(), null_vector_map_);
  EXPECT_EQ(mp_1_eigen.GetQPerSter(), null_vector_map_);
  EXPECT_EQ(mp_1_eigen.GetNuSigF(), correct_nu_sig_f_map_1_eigen);
  EXPECT_EQ(mp_1_eigen.GetSigS(), correct_sigma_s_map_1_eigen);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1_eigen.GetSigSPerSter(), correct_sigma_s_per_ster_map_1_eigen);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1_eigen.GetChiNuSigF(), correct_chi_nu_sig_f_map_1_eigen);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_1_eigen.GetChiNuSigFPerSter(), correct_chi_nu_sig_f_per_ster_map_1_eigen);}

  // four material eigen problem
  map[11] = test_materials_.at(11);
  const std::vector<double> material_11_fake_chi({0.4, 0.3, 0.1, 0.2, 0, 0, 0});
  GetPropertyPtr(map.at(11), Material::CHI)->clear_value();
  for (const double& val : material_11_fake_chi) {
    GetPropertyPtr(map.at(11), Material::CHI)->add_value(val);
  }
  correct_fissile_id_map[11] = true;
  correct_sig_t_map[11] = test_sigma_t_.at(11);
  correct_inv_sig_t_map[11] =
   {8.99847026005579, 3.4669255304396063, 2.203516812833282,
    2.15526531316005, 1.453594011192674, 1.0109281331190165, 0.5949193884228687};
  correct_q_map[11] = {0, 0, 0, 0, 0, 0, 0};
  correct_q_per_ster_map[11] = {0, 0, 0, 0, 0, 0, 0};
  sigma_s_values[11] = test_sigma_s_.at(11);
  sigma_s_per_ster_values[11] = {
    0.009739486742508536, 0.005342035664879467, 2.288648081661455e-05, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.03421194656703382, 0.004110653870177472, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.059079906425142464, 0.003888712302035825, 0.0005468165957279794, 9.895458586738594e-05, 6.109401222997041e-05,
    0.0, 0.0, 0.0, 0.04355036285295078, 0.015849445007806398, 0.002175886804480847, 0.0009729141671207563,
    0.0, 0.0, 0.0, 0.0003754862994895541, 0.053277117200011964, 0.01666511409115236, 0.005174763819689885,
    0.0, 0.0, 0.0, 0.0, 0.011292838987085432, 0.06440523082099728, 0.020248487634866384,
    0.0, 0.0, 0.0, 0.0, 0.004728493359260211, 0.032808995743678765, 0.10521733287805202};

  correct_nu_sig_f_map[11] = test_nu_sig_f_.at(11);
  chi_nu_sig_f_values[11] = {
  0.0045832, 0.0034373999999999997, 0.0011458, 0.0022916, 0.0, 0.0, 0.0,
  0.00042160000000000006,  0.0003162,  0.00010540000000000001,  0.00021080000000000003,  0.0,  0.0,  0.0,
  0.004920000000000001,  0.0036899999999999997,  0.0012300000000000002,  0.0024600000000000004,  0.0,  0.0,  0.0,
  0.0090404, 0.0067802999999999995, 0.0022601, 0.0045202, 0.0, 0.0, 0.0,
  0.0383972, 0.028797899999999998, 0.0095993, 0.0191986, 0.0, 0.0, 0.0,
  0.063544, 0.047658, 0.015886, 0.031772, 0.0, 0.0, 0.0,
  0.118224, 0.088668, 0.029556, 0.059112, 0.0, 0.0, 0.0};

  chi_nu_sig_f_per_ster_values[11] = {
  0.00036471946758938735, 0.0002735396006920405, 9.117986689734684e-05, 0.00018235973379469368, 0.0, 0.0, 0.0,
  3.3549862003771544e-05, 2.516239650282865e-05, 8.387465500942886e-06, 1.6774931001885772e-05, 0.0, 0.0, 0.0,
  0.0003915211600060626, 0.0002936408700045469, 9.788029000151564e-05, 0.0001957605800030313, 0.0, 0.0, 0.0,
  0.0007194121737639854, 0.000539559130322989, 0.00017985304344099635, 0.0003597060868819927, 0.0, 0.0, 0.0,
  0.003055552090444062, 0.002291664067833046, 0.0007638880226110155, 0.001527776045222031, 0.0, 0.0, 0.0,
  0.005056670851915699, 0.003792503138936774, 0.0012641677129789247, 0.0025283354259578493, 0.0, 0.0, 0.0,
  0.009407966996048117, 0.007055975247036088, 0.0023519917490120294, 0.004703983498024059, 0.0, 0.0, 0.0};

  map[2] = test_materials_.at(2);

  correct_fissile_id_map[2] = false;
  correct_sig_t_map[2] = test_sigma_t_.at(2);
  correct_inv_sig_t_map[2] =
   {13.265414411546216, 4.0205853972338375, 2.3717477409102767,
    1.880300095895305, 1.1007275809309953, 0.7572889057175313, 0.4317230065190174};
  correct_q_map[2] = {0, 0, 0, 0, 0, 0, 0};
  correct_q_per_ster_map[2] = {0, 0, 0, 0, 0, 0, 0};
  sigma_s_values[2] = test_sigma_s_.at(2);

  sigma_s_per_ster_values[2] = {
    0.006582330136394607, 0.006522408300320508, 4.109539785575829e-05, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.03751520741090611, 0.007936261237277361, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.07603786561158392, 0.008764662716070678, 0.0012534247543202217, 0.00022825206163524172, 0.00014214127967537174,
    0.0, 0.0, 0.0, 0.056272413229001436, 0.02817758689970461, 0.003975531323492453, 0.0017766466297348276,
    0.0, 0.0, 0.0, 0.00018987980485578573, 0.0709854919431317, 0.03409496768386128, 0.010807416410655152,
    0.0, 0.0, 0.0, 0.0, 0.01724682540815324, 0.09500754327870692, 0.034432376163216094,
    0.0, 0.0, 0.0, 0.0, 0.007817452072259262, 0.05136327900933192, 0.1610090981789159};

  correct_sigma_s_map[2] = dealii::FullMatrix<double>(7, 7, sigma_s_values[2].data());
  correct_sigma_s_per_ster_map[2] = dealii::FullMatrix<double>(7, 7, sigma_s_per_ster_values[2].data());
  correct_sigma_s_map[11] = dealii::FullMatrix<double>(7, 7, sigma_s_values[11].data());
  correct_sigma_s_per_ster_map[11] = dealii::FullMatrix<double>(7, 7, sigma_s_per_ster_values[11].data());
  correct_chi_nu_sig_f_map[11] = dealii::FullMatrix<double>(7, 7, chi_nu_sig_f_values[11].data());
  correct_chi_nu_sig_f_per_ster_map[11] = dealii::FullMatrix<double>(7, 7, chi_nu_sig_f_per_ster_values[11].data());

  MaterialProtobuf mp_4_eigen(map, true, false, 7, 4);

  EXPECT_EQ(mp_4_eigen.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_4_eigen.GetSigT(), correct_sig_t_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_4_eigen.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_4_eigen.GetQ(), null_vector_map_);
  EXPECT_EQ(mp_4_eigen.GetQPerSter(), null_vector_map_);
  EXPECT_EQ(mp_4_eigen.GetNuSigF(), correct_nu_sig_f_map);
  EXPECT_EQ(mp_4_eigen.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_4_eigen.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_4_eigen.GetChiNuSigF(), correct_chi_nu_sig_f_map);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_4_eigen.GetChiNuSigFPerSter(), correct_chi_nu_sig_f_per_ster_map);}

  // four material fixed source problem
  correct_fissile_id_map[10] = false;
  correct_fissile_id_map[11] = false;
  map[10].set_is_fissionable(false);
  map[11].set_is_fissionable(false);
  
  MaterialProtobuf mp_4_q(map, false, false, 7, 4);

  EXPECT_EQ(mp_4_q.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_4_q.GetSigT(), correct_sig_t_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_4_q.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_4_q.GetQ(), correct_q_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_4_q.GetQPerSter(), correct_q_per_ster_map);};
  EXPECT_EQ(mp_4_q.GetNuSigF(), null_vector_map_);
  EXPECT_EQ(mp_4_q.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_4_q.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  EXPECT_EQ(mp_4_q.GetChiNuSigF(), null_matrix_map_);
  EXPECT_EQ(mp_4_q.GetChiNuSigFPerSter(), null_matrix_map_);

  // one group problem
  const std::unordered_set<int> material_ids = {1, 2, 10, 11};
  for (int key : material_ids) {
    correct_sig_t_map[key].resize(1);
    correct_inv_sig_t_map[key].resize(1);
    correct_q_map[key].resize(1);
    correct_q_per_ster_map[key].resize(1);
    correct_sigma_s_map[key] = dealii::FullMatrix<double>(1, 1);
    correct_sigma_s_map[key].set(0, 0, sigma_s_values[key][0]);
    correct_sigma_s_per_ster_map[key] = dealii::FullMatrix<double>(1, 1);
    correct_sigma_s_per_ster_map[key].set(0, 0, sigma_s_per_ster_values[key][0]);
    if (correct_nu_sig_f_map.count(key) > 0) {
      correct_nu_sig_f_map[key].resize(1);
      correct_chi_nu_sig_f_map[key] = dealii::FullMatrix<double>(1, 1);
      correct_chi_nu_sig_f_map[key].set(0, 0, correct_nu_sig_f_map[key][0]);
      correct_chi_nu_sig_f_per_ster_map[key] = dealii::FullMatrix<double>(1, 1);
      correct_chi_nu_sig_f_per_ster_map[key].set(0, 0, correct_nu_sig_f_map[key][0] / (4*3.141592653589793));
    }
  
    map[key].set_number_of_groups(1);
    GetPropertyPtr(map[key], Material::ENERGY_GROUPS)->clear_value();
    GetPropertyPtr(map[key], Material::ENERGY_GROUPS)->add_value(20e6);
    GetPropertyPtr(map[key], Material::ENERGY_GROUPS)->add_value(1e6);
    GetPropertyPtr(map[key], Material::SIGMA_T)->clear_value();
    GetPropertyPtr(map[key], Material::SIGMA_T)->add_value(correct_sig_t_map[key][0]);
    GetPropertyPtr(map[key], Material::SIGMA_S)->clear_value();
    GetPropertyPtr(map[key], Material::SIGMA_S)->add_value(correct_sigma_s_map[key](0, 0));
    if (correct_nu_sig_f_map.count(key) > 0) {
      GetPropertyPtr(map[key], Material::NU_SIG_F)->clear_value();
      GetPropertyPtr(map[key], Material::NU_SIG_F)->add_value(correct_nu_sig_f_map[key][0]);
    }
    if (GetPropertyPtr(map[key], Material::CHI) != nullptr) {
      GetPropertyPtr(map[key], Material::CHI)->clear_value();
      GetPropertyPtr(map[key], Material::CHI)->add_value(1);
    }
    if (GetPropertyPtr(map[key], Material::Q) != nullptr) {
      GetPropertyPtr(map[key], Material::Q)->clear_value();
      GetPropertyPtr(map[key], Material::Q)->add_value(correct_q_map[key][0]);
    }
  }

  MaterialProtobuf mp_one_group(map, false, false, 1, 4);

  EXPECT_EQ(mp_one_group.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_one_group.GetSigT(), correct_sig_t_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_one_group.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_one_group.GetQ(), correct_q_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_one_group.GetQPerSter(), correct_q_per_ster_map);};
  EXPECT_EQ(mp_one_group.GetNuSigF(), null_vector_map_);
  EXPECT_EQ(mp_one_group.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_one_group.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  EXPECT_EQ(mp_one_group.GetChiNuSigF(), null_matrix_map_);
  EXPECT_EQ(mp_one_group.GetChiNuSigFPerSter(), null_matrix_map_);

  // one group eigen problem
  correct_fissile_id_map[10] = true;
  correct_fissile_id_map[11] = true;
  map[10].set_is_fissionable(true);
  map[11].set_is_fissionable(true);
  MaterialProtobuf mp_one_group_eigen(map, true, false, 1, 4);

  EXPECT_EQ(mp_one_group_eigen.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_one_group_eigen.GetSigT(), correct_sig_t_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_one_group_eigen.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_one_group_eigen.GetQ(), null_vector_map_);
  EXPECT_EQ(mp_one_group_eigen.GetQPerSter(), null_vector_map_);
  EXPECT_EQ(mp_one_group_eigen.GetNuSigF(), correct_nu_sig_f_map);
  EXPECT_EQ(mp_one_group_eigen.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_one_group_eigen.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_one_group_eigen.GetChiNuSigF(), correct_chi_nu_sig_f_map);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_one_group_eigen.GetChiNuSigFPerSter(), correct_chi_nu_sig_f_per_ster_map);}
}

TEST_F(MaterialProtobufTest, FailedToFindMaterialFile) {
  dealii::ParameterHandler prm;
  prm.declare_entry("do eigenvalue calculations", "false", dealii::Patterns::Bool());
  prm.declare_entry("do nda", "false", dealii::Patterns::Bool());
  prm.declare_entry("number of groups", "7", dealii::Patterns::Integer());
  prm.declare_entry("number of materials", "1", dealii::Patterns::Integer());
  prm.enter_subsection("material ID map");
  prm.declare_entry("material id file name map", "3 : src/material/tests/data/broken/nonexistent_file_70844.material",
    dealii::Patterns::Map(dealii::Patterns::Integer(), dealii::Patterns::FileName(dealii::Patterns::FileName::input)));
  prm.leave_subsection();
  prm.enter_subsection("fissile material IDs");
  prm.declare_entry("fissile material ids", "", 
    dealii::Patterns::List(dealii::Patterns::Integer()));
  prm.leave_subsection();

  EXPECT_THROW({
    try {
      MaterialProtobuf mp(prm);
    }
    catch (const MaterialProtobuf::FailedToFindMaterialFile& e) {
      std::string expected = "Failed to find material file \"src/material/tests/data/broken/nonexistent_file_70844.material\"";
      expected += " for material number 3 in the current working directory.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::FailedToFindMaterialFile);

  prm.set("number of materials", "2");
  prm.enter_subsection("material ID map");
  prm.set("material id file name map",
    "10 : src/material/tests/data/serialized/uo2_20.material, 6 : src/material/tests/data/broken/nonexistent_file_55314.material");
  prm.leave_subsection();

  EXPECT_THROW({
    try {
      MaterialProtobuf mp(prm);
    }
    catch (const MaterialProtobuf::FailedToFindMaterialFile& e) {
      std::string expected = "Failed to find material file \"src/material/tests/data/broken/nonexistent_file_55314.material\"";
      expected += " for material number 6 in the current working directory.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::FailedToFindMaterialFile);

}

TEST_F(MaterialProtobufTest, FailedToParseMaterialFile) {
  dealii::ParameterHandler prm;
  prm.declare_entry("do eigenvalue calculations", "false", dealii::Patterns::Bool());
  prm.declare_entry("do nda", "false", dealii::Patterns::Bool());
  prm.declare_entry("number of groups", "7", dealii::Patterns::Integer());
  prm.declare_entry("number of materials", "1", dealii::Patterns::Integer());
  prm.enter_subsection("material ID map");
  prm.declare_entry("material id file name map", "3 : ./test_data/material/broken/invalid_material.material",
    dealii::Patterns::Map(dealii::Patterns::Integer(), dealii::Patterns::FileName(dealii::Patterns::FileName::input)));
  prm.leave_subsection();
  prm.enter_subsection("fissile material IDs");
  prm.declare_entry("fissile material ids", "", 
    dealii::Patterns::List(dealii::Patterns::Integer()));
  prm.leave_subsection();

  EXPECT_THROW({
    try {
      /*
        As long as the LogSilencer object exists, non-fatal protobuf errors will be discarded.
        Otherwise, problems with parsing our intentionally broken file would be printed to stderr.
      */
      google::protobuf::LogSilencer log_silencer;
      MaterialProtobuf mp(prm);
    }
    catch (const MaterialProtobuf::FailedToParseMaterialFile& e) {
      std::string expected = "Failed to parse file \"./test_data/material/broken/invalid_material.material\"";
      expected += " for material number 3 as either a human-readable or serialized material file defined by material.proto.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::FailedToParseMaterialFile);

  prm.set("number of materials", "2");
  prm.enter_subsection("material ID map");
  prm.set("material id file name map",
    "10 : ./test_data/material/serialized/uo2_20.material, 26 : ./test_data/material/broken/invalid_material.material");
  prm.leave_subsection();

  EXPECT_THROW({
    try {
      google::protobuf::LogSilencer log_silencer;
      MaterialProtobuf mp(prm);
    }
    catch (const MaterialProtobuf::FailedToParseMaterialFile& e) {
      std::string expected = "Failed to parse file \"./test_data/material/broken/invalid_material.material\"";
      expected += " for material number 26 as either a human-readable or serialized material file defined by material.proto.";
      EXPECT_EQ(expected, GetMessage(e));
      throw;
    }
  }, MaterialProtobuf::FailedToParseMaterialFile);
}

TEST_F(MaterialProtobufTest, ConstructorFromParameterHandler) {
  dealii::ParameterHandler prm;
  prm.declare_entry("do eigenvalue calculations", "false", dealii::Patterns::Bool());
  prm.declare_entry("do nda", "false", dealii::Patterns::Bool());
  prm.declare_entry("number of groups", "7", dealii::Patterns::Integer());
  prm.declare_entry("number of materials", "4", dealii::Patterns::Integer());
  prm.enter_subsection("fissile material IDs");
  prm.declare_entry("fissile material ids", "", 
    dealii::Patterns::List(dealii::Patterns::Integer()));
  prm.leave_subsection();

  // try two in serialized format and two in human readable text format
  prm.enter_subsection("material ID map");
  const std::string path = "./test_data/material/";
  std::string entry = "1 : " + path + "serialized/control_rod.material, ";
  entry += "11 : " + path + "readable/uo2_33.material, ";
  entry += "10: " + path + "serialized/uo2_20.material, ";
  entry += "2: " + path + "readable/reflector.material";
  prm.declare_entry("material id file name map", entry,
    dealii::Patterns::Map(dealii::Patterns::Integer(), dealii::Patterns::FileName(dealii::Patterns::FileName::input)));
  prm.leave_subsection();

  std::unordered_set<int> material_ids = {1, 2, 10, 11};
  std::unordered_map<int, bool> correct_fissile_id_map = {{1, false}, {2, false}, {10, false}, {11, false}};
  std::unordered_map<int, std::vector<double>> zero_vals_map = 
  {{1, {0, 0, 0, 0, 0, 0, 0}}, {2, {0, 0, 0, 0, 0, 0, 0}}, {10, {0, 0, 0, 0, 0, 0, 0}}, {11, {0, 0, 0, 0, 0, 0, 0}}};

  std::unordered_map<int, std::vector<double>> correct_inv_sig_t_map;
  correct_inv_sig_t_map[1] =
   {8.053475074494644, 3.34213428695565, 1.7242271151956137,
    0.9450902561194594, 0.7574036203893054, 0.6134592969756456, 0.43769422681314835};
  correct_inv_sig_t_map[10] =
   {9.000090000900009, 3.464643314970724, 2.2173932325158545,
    2.1791714790036827, 1.4955954713369128, 1.0481960545900506, 0.6233248145608676};
  correct_inv_sig_t_map[11] =
   {8.99847026005579, 3.4669255304396063, 2.203516812833282,
    2.15526531316005, 1.453594011192674, 1.0109281331190165, 0.5949193884228687};
  correct_inv_sig_t_map[2] =
   {13.265414411546216, 4.0205853972338375, 2.3717477409102767,
    1.880300095895305, 1.1007275809309953, 0.7572889057175313, 0.4317230065190174};

  std::unordered_map<int, std::vector<double>> sigma_s_per_ster_values;
  sigma_s_per_ster_values[1] = {
    0.009815085340477186, 0.00442490530531242, 1.8601233973865266e-05, 0, 0, 0, 0,
    0, 0.032165213998872053, 0.00345541296946814, 0, 0, 0, 0,
    0, 0, 0.052924589001063414, 0.0033144812673602665, 0.00046515419442752796, 8.417704940130345e-05, 5.197045511722751e-05,
    0, 0, 0, 0.04337370086611878, 0.016879973264326418, 0.0023399755508085912, 0.001046364173357666,
    0, 0, 0, 0.0003438303813085761, 0.0547827229616613, 0.018194593094265476, 0.005781144152870007,
    0, 0, 0, 0, 0.011935029182461232, 0.06806101986381811, 0.02220529766018124,
    0, 0, 0, 0, 0.005154312409502577, 0.03487482690501156, 0.11221219262694081};
    sigma_s_per_ster_values[2] = {
    0.006582330136394607, 0.006522408300320508, 4.109539785575829e-05, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.03751520741090611, 0.007936261237277361, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.07603786561158392, 0.008764662716070678, 0.0012534247543202217, 0.00022825206163524172, 0.00014214127967537174,
    0.0, 0.0, 0.0, 0.056272413229001436, 0.02817758689970461, 0.003975531323492453, 0.0017766466297348276,
    0.0, 0.0, 0.0, 0.00018987980485578573, 0.0709854919431317, 0.03409496768386128, 0.010807416410655152,
    0.0, 0.0, 0.0, 0.0, 0.01724682540815324, 0.09500754327870692, 0.034432376163216094,
    0.0, 0.0, 0.0, 0.0, 0.007817452072259262, 0.05136327900933192, 0.1610090981789159};
  sigma_s_per_ster_values[10] = {
    0.009743465616085833, 0.00534410467913966173, 2.288648081661455e-05, 0, 0, 0, 0,
    0, 0.03423422825906669, 0.004111370067421386, 0, 0, 0, 0,
    0, 0, 0.05904727966180864, 0.0038879165273203653, 0.0005463709618873221, 9.887500839583999e-05, 6.104467419761192e-05,
    0, 0, 0, 0.04348829242514494, 0.01583512106292813, 0.0021727832830905555, 0.00097156135010447526,
    0, 0, 0, 0.00036002439676817644, 0.05301689886805672, 0.016649198596843173, 0.005144524380502425,
    0, 0, 0, 0, 0.011172677005051052, 0.06390150542611144, 0.020053522829578813,
    0, 0, 0, 0, 0.0046220982798032785, 0.032149298504562863, 0.10337113553818602};
  sigma_s_per_ster_values[11] = {
    0.009739486742508536, 0.005342035664879467, 2.288648081661455e-05, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.03421194656703382, 0.004110653870177472, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.059079906425142464, 0.003888712302035825, 0.0005468165957279794, 9.895458586738594e-05, 6.109401222997041e-05,
    0.0, 0.0, 0.0, 0.04355036285295078, 0.015849445007806398, 0.002175886804480847, 0.0009729141671207563,
    0.0, 0.0, 0.0, 0.0003754862994895541, 0.053277117200011964, 0.01666511409115236, 0.005174763819689885,
    0.0, 0.0, 0.0, 0.0, 0.011292838987085432, 0.06440523082099728, 0.020248487634866384,
    0.0, 0.0, 0.0, 0.0, 0.004728493359260211, 0.032808995743678765, 0.10521733287805202};

  std::unordered_map<int, std::vector<double>> correct_nu_sig_f_map;
  correct_nu_sig_f_map[10] = test_nu_sig_f_[10];
  correct_nu_sig_f_map[11] = test_nu_sig_f_[11];

  std::unordered_map<int, std::vector<double>> chi_nu_sig_f_values;
  chi_nu_sig_f_values[10] = {
   0.006671236516878692, 0.004584026253091862, 3.737230029446337e-06, 0, 0, 0, 0,
   0.00040727190868430316, 0.00027984993738899475, 2.2815392670218845e-07, 0, 0, 0, 0,
   0.004584245730863049, 0.0031499861723884113, 2.5680967485407607e-06, 0, 0, 0, 0,
   0.008204690651853561, 0.005637713076344526, 4.596271801913441e-06, 0, 0, 0, 0,
   0.03563562541949714, 0.024486411485340687, 1.9963095162177958e-05, 0, 0, 0, 0,
   0.05847508565394124, 0.040180156573863555, 3.275777219522161e-05, 0, 0, 0, 0,
   0.11205828671037378, 0.07699893818143114, 6.277510819512313e-05, 0, 0, 0, 0};
  chi_nu_sig_f_values[11] = {
   0.006789148948432014, 0.004665047766935479, 3.8032846325069838e-06, 0.0, 0.0, 0.0, 0.0,
   0.0006245211198854376, 0.0004291290230712162, 3.4985704334633976e-07, 0.0, 0.0, 0.0, 0.0,
   0.007288054814602356, 0.0050078624134496765, 4.08277194796962e-06, 0.0, 0.0, 0.0, 0.0,
   0.013391652590636408, 0.009201845398892369, 7.502010471224502e-06, 0.0, 0.0, 0.0, 0.0,
   0.05687823136732714, 0.03908290541904673, 3.186321362613396e-05, 0.0, 0.0, 0.0, 0.0,
   0.09412848681688864, 0.06467878235777363, 5.273082533776047e-05, 0.0, 0.0, 0.0, 0.0,
   0.1751266244718595, 0.12033526950562491, 9.810602251560169e-05, 0.0, 0.0, 0.0, 0.0};
  std::unordered_map<int, std::vector<double>> chi_nu_sig_f_per_ster_values;
  chi_nu_sig_f_per_ster_values[10] = {
   0.0005308801340982011, 0.0003647852187212948, 2.973993163289271e-07, 0.0, 0.0, 0.0, 0.0,
   3.2409668724788937e-05, 2.2269750429707967e-05, 1.815591261023963e-08, 0.0, 0.0, 0.0, 0.0,
   0.00036480268420738634, 0.0002506679350033674, 2.043626459342431e-07, 0.0, 0.0, 0.0, 0.0,
   0.0006529085368911796, 0.00044863495191702367, 3.657596885342085e-07, 0.0, 0.0, 0.0, 0.0,
   0.0028357929678420836, 0.0019485667132370648, 1.588612637237262e-06, 0.0, 0.0, 0.0, 0.0,
   0.004653299464773362, 0.0031974352664683493, 2.6067806847738836e-06, 0.0, 0.0, 0.0, 0.0,
   0.008917315122182416, 0.00612738081220102, 4.995484386191197e-06, 0.0, 0.0, 0.0, 0.0};
  chi_nu_sig_f_per_ster_values[11] = {
   0.0005402633072650492, 0.0003712327059337948, 3.0265577462446457e-07, 0.0, 0.0, 0.0, 0.0,
   4.9697811647526785e-05, 3.414900262299003e-05, 2.7840738912042733e-08, 0.0, 0.0, 0.0, 0.0,
   0.0005799649746343259, 0.00039851302871231245, 3.2489666851814576e-07, 0.0, 0.0, 0.0, 0.0,
   0.0010656738529845851, 0.000732259590400567, 5.969910248112692e-07, 0.0, 0.0, 0.0, 0.0,
   0.004526225838217304, 0.003110118793917155, 2.535593975696127e-06, 0.0, 0.0, 0.0, 0.0,
   0.007490506981334066, 0.005146973962702274, 4.19618575290997e-06, 0.0, 0.0, 0.0, 0.0,
   0.013936133975847262, 0.009575976485057812, 7.807029215221395e-06, 0.0, 0.0, 0.0, 0.0};

  std::unordered_map<int, dealii::FullMatrix<double>> correct_sigma_s_map;
  std::unordered_map<int, dealii::FullMatrix<double>> correct_sigma_s_per_ster_map;
  std::unordered_map<int, dealii::FullMatrix<double>> correct_chi_nu_sig_f_map;
  std::unordered_map<int, dealii::FullMatrix<double>> correct_chi_nu_sig_f_per_ster_map;
  for (int id : material_ids) {
    correct_sigma_s_map[id] = dealii::FullMatrix<double>(7, 7, test_sigma_s_[id].data());
    correct_sigma_s_per_ster_map[id] = dealii::FullMatrix<double>(7, 7, sigma_s_per_ster_values[id].data());
    if (chi_nu_sig_f_values.count(id) > 0) {
      correct_chi_nu_sig_f_map[id] = dealii::FullMatrix<double>(7, 7, chi_nu_sig_f_values[id].data());
      correct_chi_nu_sig_f_per_ster_map[id] = dealii::FullMatrix<double>(7, 7, chi_nu_sig_f_per_ster_values[id].data());
    }
  }

  std::stringstream prm_string_before;
  prm.print_parameters(prm_string_before, dealii::ParameterHandler::Text);

  MaterialProtobuf mp(prm);

  std::stringstream prm_string_after;
  prm.print_parameters(prm_string_after, dealii::ParameterHandler::Text);
  EXPECT_EQ(prm_string_before.str(), prm_string_after.str());

  EXPECT_EQ(mp.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp.GetSigT(), test_sigma_t_);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp.GetQ(), zero_vals_map);
  EXPECT_EQ(mp.GetQPerSter(), zero_vals_map);
  EXPECT_EQ(mp.GetNuSigF(), null_vector_map_);
  EXPECT_EQ(mp.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  EXPECT_EQ(mp.GetChiNuSigF(), null_matrix_map_);
  EXPECT_EQ(mp.GetChiNuSigFPerSter(), null_matrix_map_);

  prm.set("do eigenvalue calculations", "true");
  prm.enter_subsection("fissile material IDs");
  prm.set("fissile material ids", "11, 10");
  prm.leave_subsection();
  correct_fissile_id_map[10] = true;
  correct_fissile_id_map[11] = true;
  
  prm.print_parameters(prm_string_before, dealii::ParameterHandler::Text);

  MaterialProtobuf mp_eigen(prm);

  prm.print_parameters(prm_string_after, dealii::ParameterHandler::Text);
  EXPECT_EQ(prm_string_before.str(), prm_string_after.str());
  // make sure the parameter handler path is back where it started
  EXPECT_NO_THROW(prm.set("do eigenvalue calculations", "false"));

  EXPECT_EQ(mp_eigen.GetFissileIDMap(), correct_fissile_id_map);
  EXPECT_EQ(mp_eigen.GetSigT(), test_sigma_t_);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_eigen.GetInvSigT(), correct_inv_sig_t_map);}
  EXPECT_EQ(mp_eigen.GetQ(), null_vector_map_);
  EXPECT_EQ(mp_eigen.GetQPerSter(), null_vector_map_);
  EXPECT_EQ(mp_eigen.GetNuSigF(), correct_nu_sig_f_map);
  EXPECT_EQ(mp_eigen.GetSigS(), correct_sigma_s_map);
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_eigen.GetSigSPerSter(), correct_sigma_s_per_ster_map);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_eigen.GetChiNuSigF(), correct_chi_nu_sig_f_map);}
  {SCOPED_TRACE("");
  ExpectApproxEqual(mp_eigen.GetChiNuSigFPerSter(), correct_chi_nu_sig_f_per_ster_map);}
}

TEST_F(MaterialProtobufTest, DiffusionCoefficientValuesTest) {

  auto &diffusion_material = test_materials_.at(1);

  auto diffusion_ptr = diffusion_material.add_vector_property();
  
  diffusion_ptr->set_id(Material::DIFFUSION_COEFF);
  
  const std::vector<double> diffusion_coef(test_helpers::RandomVector(7, 1e-3, 3));
  const std::vector<double> null_vector(7,0);
  
  for (const double& val : diffusion_coef) {
    diffusion_ptr->add_value(val);
  }
   
  MaterialProtobuf mp_1(test_materials_, false, false, 7, 4);

  auto diffusion_map = mp_1.GetDiffusionCoef();
  
  ASSERT_EQ(diffusion_map.size(), 4)
      << "Incorrectly sized diffusion map";
  ASSERT_EQ(mp_1.GetDiffusionCoef().at(1), diffusion_coef)
      << "Incorrect provided diffusion coefficient";
  ASSERT_EQ(mp_1.GetDiffusionCoef().at(10), null_vector)
      << "Incorrect default diffusion coefficient";
}

TEST_F(MaterialProtobufTest, DiffusionCoefficientInvalidTest) {

  auto &diffusion_material = test_materials_.at(1);

  auto diffusion_ptr = diffusion_material.add_vector_property();
  
  diffusion_ptr->set_id(Material::DIFFUSION_COEFF);
  
  const std::vector<double> diffusion_coef(test_helpers::RandomVector(5, 1e-3, 3));
  for (const double& val : diffusion_coef) {
    diffusion_ptr->add_value(val);
  }
  
  ASSERT_THROW({
      MaterialProtobuf mp_1(test_materials_, false, false, 7, 4);
    }, MaterialProtobuf::WrongNumberOfValues);
}

} // namespace