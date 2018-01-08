#include "../aqdata/aq_base.h"

#ifdef TEST
#include ""


int main ()
{
#ifdef TEST
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TEST;
#else
  // TODO: modifications needed along with the coding progress.
#endif
}
