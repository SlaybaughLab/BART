#ifndef BART_SRC_FORMULATION_STAMPER_I_H_
#define BART_SRC_FORMULATION_STAMPER_I_H_

namespace bart {

namespace formulation {

/*! \brief Interface for a system matrix stamper.
 * See the default implementation for detailed description for the purpose of
 * this class.
 *
 * \tparam dim spatial dimension of the cells in the mesh.
 *
 * \author J.S. Rehak
 */
template <int dim>
class StamperI {
 public:
  virtual ~StamperI() = default;
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_STAMPER_I_H_