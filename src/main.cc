#ifdef TEST
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <unistd.h>
#include <cstdlib>
//bool TEST_REPORT = false;
#endif

//#include "aqdata/aq_base.h"

int main(int argc, char* argv[]) {
#ifdef TEST
  // Parse optional arguments
  // int c;
  // while ((c = getopt (argc, argv, "r")) != -1)
  //   switch(c) {
  //     case 'r':
  //       TEST_REPORT = true;
  //   }
        
  // // Testing
  ::testing::InitGoogleMock(&argc, argv);
  return RUN_ALL_TESTS();
#else  
  return 0;
#endif
}
