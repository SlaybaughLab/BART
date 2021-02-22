#include "cross_sections.hpp"

namespace bart::data::cross_sections {

CrossSections::CrossSections(material::MaterialI &materials)
    : diffusion_coef_(materials.GetDiffusionCoef()),
      sigma_t_(materials.GetSigT()),
      inverse_sigma_t_(materials.GetInvSigT()),
      sigma_s_(materials.GetSigS()),
      sigma_s_per_ster_(materials.GetSigSPerSter()),
      q_(materials.GetQ()),
      q_per_ster_(materials.GetQPerSter()),
      is_material_fissile_(materials.GetFissileIDMap()),
      nu_sigma_f_(materials.GetNuSigF()),
      fiss_transfer_(materials.GetChiNuSigF()),
      fiss_transfer_per_ster_(materials.GetChiNuSigFPerSter())
{}

} // namespace bart::data::cross_sections