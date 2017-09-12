#ifndef __problem_definition_h__
#define __problem_definition_h__

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/quadrature_lib.h>

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>

using namespace dealii;

class ProblemDefinition
{
public:
  ProblemDefinition (ParameterHandler &prm);
  ~ProblemDefinition ();

  static void declare_parameters (ParameterHandler &prm);
  
private:
  const unsigned int nmat;
  const unsigned int ngrp;
  const unsigned int z_levels;
  const unsigned int y_levels;
};


#endif  // define  __properties_h__
