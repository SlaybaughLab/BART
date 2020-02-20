#ifndef BART_SRC_FORMULATION_STAMPER_H_
#define BART_SRC_FORMULATION_STAMPER_H_

#include <memory>

#include "domain/definition_i.h"
#include "formulation/stamper_i.h"

namespace bart {

namespace formulation {

/*! \brief Stamps a system matrix or vector with a provided function over all cells.
 *  The system matrix and vectors represent the degrees of freedom for all the
 *  cells in the triangulation. This stamper class takes a cell-based function
 *  and fills the system matrix with the results of the function over all the
 *  cells in the triangulation. Another class provides a functional that will
 *  fill a matrix or vector for a given cell, this class then iterates over all
 *  cells and stamps the results onto the system matrix.
 *
 *  \author J.S. Rehak
 *  \tparam dim spatial dimension of the cells in the mesh
 */
template <int dim>
class Stamper : public StamperI<dim> {
 public:
  /*! \brief Constructor.
   * Takes a domain definition dependency that provides the domain of cells to
   * iterate over. The matrices and vectors passed in this classes functions
   * should be seperately initialized using this domain.
   */
  explicit Stamper(std::shared_ptr<domain::DefinitionI<dim>>);
  virtual ~Stamper() = default;

  void StampMatrix(
      system::MPISparseMatrix& to_stamp,
      std::function<void(formulation::FullMatrix&,
                         const domain::CellPtr<dim>&)> stamp_function)
  override;

  void StampVector(
      system::MPIVector& to_stamp,
      std::function<void(formulation::Vector&,
                         const domain::CellPtr<dim>&)> stamp_function)
  override;

  /*! \brief Access domain definition dependency */
  domain::DefinitionI<dim>* domain_ptr() const { return domain_ptr_.get(); }
 private:
  std::shared_ptr<domain::DefinitionI<dim>> domain_ptr_;
};

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_STAMPER_H_
