#ifndef BART_SOLVER_LINEAR_I_H_
#define BART_SOLVER_LINEAR_I_H_

namespace bart {

namespace solver {

/*! \brief Linear solver class.
 *
 * This class is all solvers that solve an equation in the form \f$Ax = b\f$.
 *
 */
class LinearI {
 public:
  virtual ~LinearI() = default;
};

} // namespace solver

} // namespace bart

#endif // BART_SOLVER_LINEAR_I_H_