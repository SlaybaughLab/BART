#include "../../src/aqdata/aq_base.h"
#include "../test_utilities.h"

void Test (int dummy_num) {
  dealii::deallog << "See test " << dummy_num << std::endl;
}

int main () {
  testing::init_log ();

  dealii::deallog.push("test test");
  Test(31415926);
  dealii::deallog.pop();
}
