#ifndef BART_SRC_FORMULATION_STAMPER_I_HPP_
#define BART_SRC_FORMULATION_STAMPER_I_HPP_

#include <functional>

#include "domain/domain_types.hpp"
#include "formulation/formulation_types.h"
#include "system/system_types.h"
#include "utility/has_description.h"

namespace bart::formulation {

/*! \brief Interface for a system matrix stamper.
 *
 *  The system matrix and vectors represent the degrees of freedom for all the
 *  cells in the triangulation. This stamper class takes a cell- or cell-and-face-based function and fills the system
 *  matrix with the results of the function over all the cells in the triangulation. Another class provides a
 *  functional that will fill a matrix or vector for a given cell, this class then iterates over all
 *  cells and stamps the results onto the system matrix.
 *
 * \tparam dim spatial dimension of the cells in the mesh.
 *
 * \author J.S. Rehak
 */
template <int dim>
class StamperI : public utility::HasDescription {
 public:
  //! A function that takes a matrix and a cell and fills the provided matrix with the cell values.
  using CellMatrixStampFunction = std::function<void(formulation::FullMatrix&, const domain::CellPtr<dim>&)>;
  //! A function that takes a vector and a cell and fills the provided vector with the cell values.
  using CellVectorStampFunction = std::function<void(formulation::Vector&, const domain::CellPtr<dim>&)>;
  //! A function that takes a matrix, a cell, and a face and fills the provided matrix with the cell face values.
  using FaceMatrixStampFunction = std::function<void(formulation::FullMatrix&, const domain::FaceIndex, const domain::CellPtr<dim>&)>;
  //! A function that takes a vector, a cell, and a face and fills the provided vector with the cell face values.
  using FaceVectorStampFunction = std::function<void(formulation::Vector&, const domain::FaceIndex, const domain::CellPtr<dim>&)>;

  virtual ~StamperI() = default;
  /*! \brief Stamp a system matrix with all cell values.
   *
   * Iterates over all cells in the triangulation, evaluates a function on that cell and the stamps the system matrix
   * with the cell values in their global DOF indices.
   *
   * @param to_stamp sparse system matrix to stamp
   * @param function function to evaluate on all cells
   */
  virtual auto StampMatrix(system::MPISparseMatrix& to_stamp, CellMatrixStampFunction function) -> void = 0;
  /*! \brief Stamp a system vector with all cell values.
   *
   * Iterates over all cells in the triangulation, evaluates a function on that cell and the stamps the system vector
   * with the cell values in their global DOF indices.
   *
   * @param to_stamp sparse system vector to stamp
   * @param function function to evaluate on all cells
   */
  virtual auto StampVector(system::MPIVector& to_stamp, CellVectorStampFunction function) -> void = 0;
  /*! \brief Stamp a system matrix with all cell face values.
   *
   * Iterates over all cells in the triangulation, evaluates a function on that cell and face and the stamps the system
   * matrix with the cell values in their global DOF indices. Note that this function will not evaluate if the cell and
   * face is on a boundary, that needs to be handled by the passed function.
   *
   * @param to_stamp sparse system matrix to stamp
   * @param function function to evaluate on all cells
   */
  virtual auto StampBoundaryMatrix(system::MPISparseMatrix& to_stamp, FaceMatrixStampFunction function) -> void = 0;
  /*! \brief Stamp a system vector with all cell face values.
   *
   * Iterates over all cells in the triangulation, evaluates a function on that cell and face and the stamps the system
   * vector with the cell values in their global DOF indices. Note that this function will not evaluate if the cell and
   * face is on a boundary, that needs to be handled by the passed function.
   *
   * @param to_stamp sparse system vector to stamp
   * @param function function to evaluate on all cells
   */
  virtual auto StampBoundaryVector(system::MPIVector& to_stamp, FaceVectorStampFunction function) -> void = 0;
};

} // namespace bart::formulation

#endif // BART_SRC_FORMULATION_STAMPER_I_HPP_