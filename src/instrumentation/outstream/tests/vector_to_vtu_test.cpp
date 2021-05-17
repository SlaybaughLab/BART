#include "instrumentation/outstream/vector_to_vtu.hpp"

#include <filesystem>

#include "domain/tests/domain_mock.hpp"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class InstrumentationOutstreamVectorToVtuTest : public ::testing::Test,
                                                public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using Domain = domain::DomainMock<dim>;

  std::shared_ptr<Domain> domain_ptr_{ nullptr };
  const std::string data_name{ "data_vector" };
  const std::string directory{ "test_output_directory" };
  const std::string filename_base{ "filename_base" };

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto InstrumentationOutstreamVectorToVtuTest<DimensionWrapper>::SetUp() -> void {
  domain_ptr_ = std::make_shared<Domain>();
  this->SetUpDealii();
}

TYPED_TEST_SUITE(InstrumentationOutstreamVectorToVtuTest, bart::testing::AllDimensions);

TYPED_TEST(InstrumentationOutstreamVectorToVtuTest, Constructor) {
  constexpr int dim = this->dim;
  instrumentation::outstream::VectorToVTU<dim> test_outstream(this->domain_ptr_,
                                                              this->data_name,
                                                              this->directory,
                                                              this->filename_base);
  EXPECT_NE(test_outstream.definition_ptr(), nullptr);
  EXPECT_EQ(test_outstream.data_name(), this->data_name);
  EXPECT_EQ(test_outstream.directory(), this->directory);
  EXPECT_EQ(test_outstream.filename_base(), this->filename_base);
}

TYPED_TEST(InstrumentationOutstreamVectorToVtuTest, Output) {
  constexpr int dim = this->dim;
  namespace fs = std::filesystem;

  instrumentation::outstream::VectorToVTU<dim> test_outstream(this->domain_ptr_, this->data_name, this->directory,
                                                              this->filename_base);

  dealii::Vector<double> output_vector(this->dof_handler_.n_dofs());
  EXPECT_CALL(*this->domain_ptr_, dof_handler()).WillOnce(::testing::ReturnRef(this->dof_handler_));
  test_outstream.Output(output_vector);
  std::string filename{ this->directory + "/" + this->filename_base + "_0000"};
  EXPECT_TRUE(fs::exists(filename + ".0.vtu"));
  EXPECT_TRUE(fs::exists(filename + ".pvtu"));
  fs::remove_all(this->directory);
}

} // namespace
