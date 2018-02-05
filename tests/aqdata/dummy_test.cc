#include "../../src/aqdata/aq_base.h"
#include "../test_utilities.h"

void test(int dummy_num)
{
  dealii::deallog << "See test " << dummy_num << std::endl;
}

int main()
{
  testing::init_log ();

  dealii::deallog.push("test test");
  test(31415926);
  dealii::deallog.pop();
}
