#ifndef BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_
#define BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_

#include "output_i.h"

namespace bart {

namespace results {

class OutputDealiiVtu : public OutputI {
 public:
  void Output(system::System &to_output) const override {};
};

} // namespace result

} //namespace bart

#endif //BART_SRC_RESULTS_OUTPUT_DEALII_VTU_H_
