#include "formulation/equation/diffusion.h"

namespace bart {

namespace formulation {

namespace equation {

template<int dim>
void Diffusion<dim>::FillCellBilinearTerm(Matrix &to_fill,
                                          const CellPtr &cell_ptr,
                                          const GroupNumber group) {
  int material_id = cell_ptr->material_id();
}

template class Diffusion<1>;
template class Diffusion<2>;
template class Diffusion<3>;

} // namespace equation

} // namespace formulation

} // namespace bart