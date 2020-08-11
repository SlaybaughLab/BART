#include "instrumentation/outstream/to_ostream.h"

#include <memory>
#include <sstream>

#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace instrumentation = bart::instrumentation;

class InstrumentationOutstreamToOstreamTest : public ::testing::Test {
 public:
  using OutstreamType = instrumentation::outstream::ToOstream;
  using OstreamType = std::ostringstream;
  std::unique_ptr<OutstreamType> test_outstream_ptr_;
  OstreamType* out_stringstream_obs_ptr_;
  void SetUp() override;
};

void InstrumentationOutstreamToOstreamTest::SetUp() {
  auto ostream_ptr = std::make_unique<OstreamType>();
  out_stringstream_obs_ptr_ = ostream_ptr.get();
  test_outstream_ptr_ = std::make_unique<OutstreamType>(
      std::move(ostream_ptr));
}

TEST_F(InstrumentationOutstreamToOstreamTest, Constructor) {
  EXPECT_NO_THROW({
    OutstreamType test_outstream(std::make_unique<OstreamType>());
  });
}

TEST_F(InstrumentationOutstreamToOstreamTest, OstreamGetter) {
  EXPECT_NE(test_outstream_ptr_->ostream_ptr(), nullptr);
}

TEST_F(InstrumentationOutstreamToOstreamTest, OutstreamString) {
  std::string test_string{"this is a test string\n"},
      second_string{"so is this!"};
  test_outstream_ptr_->Output(test_string).Output(second_string);
  EXPECT_EQ(test_string + second_string, out_stringstream_obs_ptr_->str());
}

} // namespace
