#ifndef BART_SRC_QUADRATURE_FACTORY_QUADRATURE_FACTORIES_H_
#define BART_SRC_QUADRATURE_FACTORY_QUADRATURE_FACTORIES_H_

#include <memory>

#include "quadrature/ordinate_i.h"
#include "quadrature/quadrature_generator_i.h"
#include "quadrature/quadrature_point_i.h"
#include "quadrature/quadrature_set_i.h"
#include "quadrature/quadrature_types.h"

namespace bart {

namespace quadrature {

namespace factory {

/*! \brief Factory for objects of type quadrature::OrdinateI.
 *
 * @tparam dim spatial dimension of ordinate.
 * @param type implementation type for ordiante to be created.
 * @return shared pointer to default-initialized ordinate object.
 */
template <int dim>
std::shared_ptr<OrdinateI<dim>> MakeOrdinatePtr(
    const OrdinateType type = OrdinateType::kDefault);

/*! \brief Factory for objects of type quadrature::QuadraturePointI.
 *
 * @tparam dim spatial dimension of quadrature point.
 * @param impl implementation type for quadrature point to be created.
 * @return shared pointer to default-initialized quadrature point object.
 */
template <int dim>
std::shared_ptr<QuadraturePointI<dim>> MakeQuadraturePointPtr(
    const QuadraturePointImpl impl = QuadraturePointImpl::kDefault);

/*! \brief Factory for objects of type QuadratureGeneratorI.
 *
 * @tparam dim spatial dimension of the quadrature generator.
 * @param order quadrature order.
 * @param type type of quadrature generator to be created.
 * @return shared pointer to the quadrature generator.
 */
template <int dim>
std::shared_ptr<QuadratureGeneratorI<dim>> MakeAngularQuadratureGeneratorPtr(
    const Order order, const AngularQuadratureSetType type);

/*! \brief Factory for objects of type QuadratureSetI.
 *
 * @tparam dim spatial dimension of the quadrature set.
 * @param impl implementation type for quadrature set to be created.
 * @return shared pointer to default-initialized quadrature set.
 */
template <int dim>
std::shared_ptr<QuadratureSetI<dim>> MakeQuadratureSetPtr(
    const QuadratureSetImpl impl = QuadratureSetImpl::kDefault);

/*! \brief Function to fill a quadrature set.
 *
 * The quadrature set to fill must be empty or this will throw an error.
 * Reflections across the origin will be added (unless the point is it's own
 * reflection, i.e. the origin). This is intended to be used with the
 * utility function quadrature::utility::GenerateAllPositiveX for any set not in
 * one dimension. An example is shown here, a similar procedure is used by the
 * framework builder.
 *
 * \code{.cpp}
 * std::shared_ptr<quadrature::QuadratureGeneratorI<dim>> quadrature_generator_ptr =
 *     quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(quadrature::Order(order_value),
 *                                                                 quadrature::AngularQuadratureSetType::kLevelSymmetricGaussian);
 * auto quadrature_points = quadrature::utility::GenerateAllPositiveX<dim>(quadrature_generator_ptr->GenerateSet());
 * auto quadrature_set_ptr = quadrature::factory::MakeQuadratureSetPtr<dim>();
 *
 * quadrature::factory::FillQuadratureSet<dim>(quadrature_set_ptr.get(),
                                               quadrature_points);
 * \endcode
 *
 * @tparam dim spatial dimension of quadrature set to fill.
 * @param to_fill quadrature set to fill.
 * @param point_vector vector containing points to fill the quadrature set.
 */
template <int dim>
void FillQuadratureSet(
    QuadratureSetI<dim>* to_fill,
    const std::vector<std::pair<quadrature::CartesianPosition<dim>,
                                quadrature::Weight>>& point_vector);


} // namespace factory

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_FACTORY_QUADRATURE_FACTORIES_H_
