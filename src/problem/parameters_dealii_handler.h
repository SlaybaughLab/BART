#ifndef BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_
#define BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_

#include <string>

#include <deal.II/base/parameter_handler.h>

#include "parameters_i.h"

namespace bart {

namespace problem {

/*!
 * \brief Problem parameters derived using a dealii ParameterHandler object.
 * 
 * The ParameterHandler can be given directly or a file that can be parsed into
 * a ParameterHandler object can be passed.
 */

class ParametersDealiiHandler : ParametersI {
 public:
  ParametersDealiiHandler();
  ~ParametersDealiiHandler() = default;

  /*! \brief Parses a ParameterHandler object for problem parameters */
  void Parse(const dealii::ParameterHandler &handler);
  
  /*! \brief Set up a ParameterHandler object for problem parameters.
   * This setup includes default values.
   */
  void SetUp(dealii::ParameterHandler &handler);
  
  // Basic Parameters
  /*! Get problem spatial dimension */
  std::vector<int>    NCells() const override { return n_cells_; }
  
  /*! Get problem output filename base */
  std::string         OutputFilenameBase() const override {
    return output_filename_base_; }
  
  /*! Get number of spatial dimensions */
  int                 SpatialDimension() const override {
    return spatial_dimension_; }
  
  /*! Get maximum x, y, z size */
  std::vector<double> SpatialMax() const override { return spatial_max; }


 private:
  std::vector<int>    n_cells_;
  std::string         output_filename_base_;
  int                 spatial_dimension_;
  std::vector<double> spatial_max;

  /*! Parses a ParameterHandler entry of type dealii::Patterns::List with doubles
   * into a vector.
   */
  std::vector<double> ParseDealiiList(std::string to_parse);
  std::vector<int>    ParseDealiiIntList(std::string to_parse);
  
  // Key-words for input file
  const std::string kNCells_ = "number of cells for x, y, z directions";
  const std::string kOutputFilenameBase_ = "output file name base";
  const std::string kSpatialDimension_ = "problem dimension";
  const std::string kSpatialMax_ = "x, y, z max values of boundary locations";
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_
