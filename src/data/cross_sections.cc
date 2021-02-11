#include "cross_sections.h"

namespace bart {

namespace data {

CrossSections::CrossSections(material::MaterialI &materials)
    : diffusion_coef(materials.GetDiffusionCoef()),
      sigma_t(materials.GetSigT()),
      inverse_sigma_t(materials.GetInvSigT()),
      sigma_s(materials.GetSigS()),
      sigma_s_per_ster(materials.GetSigSPerSter()),
      q(materials.GetQ()),
      q_per_ster(materials.GetQPerSter()),
      is_material_fissile(materials.GetFissileIDMap()),
      nu_sigma_f(materials.GetNuSigF()),
      fiss_transfer(materials.GetChiNuSigF()),
      fiss_transfer_per_ster(materials.GetChiNuSigFPerSter())
{}

} // namespace data

} // namespace bart 
