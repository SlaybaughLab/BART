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

const unsigned int nmat = 50;
const unsigned int ngrp = 30;
const unsigned int z_levels = 30;
const unsigned int y_levels = 100;

class ProblemDefinition
{
public:
  ProblemDefinition ();
  ~ProblemDefinition ();

  static void declare_parameters (ParameterHandler &prm);
};


#endif  // define  __properties_h__
