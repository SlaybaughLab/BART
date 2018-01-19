#include "aqdata/aq_base.h"

#ifdef TEST
#include "gtest/gtest.h"
#endif

int main(int argc, char* argv[]) {
#ifdef TEST
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
#else  
  return 0;
#endif
}
