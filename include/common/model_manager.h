#ifndef __MODEL_MANAGER__H__
#define __MODEL_MANAGER__H__

#include <deal.II/base/parameter_handler.h>

#include <string>
#include <map>

using namespace dealii;

class ModelManager
{
public:
  ModelManager (ParameterHandler &prm);
  ~ModelManager ();

  void build_and_run_model (ParameterHandler &prm);

private:
  unsigned int dim;
  std::string transport_model_name;
  std::map<std::string, unsigned int> method_index;
};

#endif //__MODEL_MANAGER__H__
