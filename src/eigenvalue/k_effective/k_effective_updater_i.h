#ifndef BART_SRC_EIGENVALUE_K_EFFECTIVE_K_EFFECTIVE_UPDATER_I_H_
#define BART_SRC_EIGENVALUE_K_EFFECTIVE_K_EFFECTIVE_UPDATER_I_H_

namespace bart {

namespace eigenvalue {

namespace k_effective {

class K_EffectiveUpdaterI {
 public:
  virtual ~K_EffectiveUpdaterI() = default;

  virtual double k_effective() const = 0;
};

} // namespace k_effective

} // namespace eigenvalue

} // namespace bart

#endif // BART_SRC_EIGENVALUE_K_EFFECTIVE_K_EFFECTIVE_UPDATER_I_H_