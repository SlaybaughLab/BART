#include "../../src/aqdata/aq_base.h"

#include <iostream>
#include <fstream>

#include <deal.II/base/logstream.h>

void test(int dummy_num)
{
  dealii::deallog << "See test " << dummy_num << std::endl;
}

int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  dealii::deallog.attach(logfile, false);
  dealii::deallog.push("test test");
  test(31415926);
  dealii::deallog.pop();
}
