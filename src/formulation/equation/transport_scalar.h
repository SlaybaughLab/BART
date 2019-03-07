#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_



#include "formulation/types.h"
#include "formulation/equation/transport.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class TransportScalar : public Transport<dim> {
 public:
  using typename Transport<dim>::CellPtr;
  using Matrix = dealii::FullMatrix<double>;
  using GroupNumber = int;

  virtual ~TransportScalar() = default;

  TransportScalar(const DiscretizationType discretization)
      : Transport<dim>(EquationType::kScalar, discretization) {};

  /*! \brief Fills a cell matrix with the integrated bilinear term for a cell.
   *
   * \param[in,out] to_fill local cell matrix to fill with the integrated term.
   * \param[in] cell_ptr the cell to use to generate the local cell matrix.
   * \param[in] group energy group number
   */
  virtual void FillCellBilinearTerm(Matrix& to_fill,
                                    const CellPtr &cell_ptr,
                                    const GroupNumber group) = 0;

  ScalarEquations scalar_equation() const { return scalar_equation_; };
 protected:
  ScalarEquations scalar_equation_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_