#include "material_cross_sections.hpp"

namespace bart::data::cross_sections {

MaterialCrossSections::MaterialCrossSections(data::material::MaterialI &materials) {
  this->diffusion_coef_ = materials.GetDiffusionCoef();
  this->sigma_t_ = materials.GetSigT();
  this->inverse_sigma_t_ = materials.GetInvSigT();
  this->sigma_s_ = materials.GetSigS();
  this->sigma_s_per_ster_ = materials.GetSigSPerSter();
  this->q_ = materials.GetQ();
  this->q_per_ster_ = materials.GetQPerSter();
  this->is_material_fissile_ = materials.GetFissileIDMap();
  this->nu_sigma_f_ = materials.GetNuSigF();
  this->fiss_transfer_ = materials.GetChiNuSigF();
  this->fiss_transfer_per_ster_ = materials.GetChiNuSigFPerSter();
}

} // namespace bart::data::cross_sections