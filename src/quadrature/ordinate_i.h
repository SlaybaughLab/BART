#ifndef BART_SRC_QUADRATURE_ORDINATE_I_H_
#define BART_SRC_QUADRATURE_ORDINATE_I_H_

namespace bart {

namespace quadrature {

/*! \brief Interface for ordinates for quadrature sets.
 *
 * @tparam dim spatial dimension.
 */
template <int dim>
class OrdinateI {
 public:
  virtual ~OrdinateI() = default;
};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_ORDINATE_I_H_
