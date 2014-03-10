/*
 * continuous_post_inters.h
 *
 *      Author: Ray
 */

#ifndef CONTINUOUS_POST_INTERS_H_
#define CONTINUOUS_POST_INTERS_H_

#include "utility/logger_stopwatch.h"

#include "core/post_operators/continuous_post.h"

#include <typeinfo>

#include "core/continuous/continuous_dynamics/ode_affine_dynamics.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/hybrid_automata/location.h"
#include "core/post_operators/sfm_post/support_function.h"

namespace hybrid_automata {

/** \brief A continuous post operator for ode_affine_dynamics.
 *
 * This version is tailored for the precise intersection algorithm;
 * it only computes the SFM for the negative invariant directions.
 * For Precise intersection, we do not intersect the
 * source invariant with the computed sfm and compute the intersection in the
 * post_discrete operator instead for precise results.
 *
 * The post_discrete operator is expected to request more directions
 * to the SFM created here.
 *
 * It is applicable under the following conditions:
 * - the dynamics are of type ode_affine_dynamics,
 * - the continuous_set class is a support_function_provider.
 */

template<class scalar_type>
class continuous_post_inters: public continuous_post {
public:
	typedef boost::shared_ptr<continuous_post_inters> ptr;
	typedef boost::shared_ptr<const continuous_post_inters> const_ptr;

	continuous_post_inters(double flowpipe_tolerance): my_flowpipe_tolerance(flowpipe_tolerance) {
	}
	;

	virtual ~continuous_post_inters() {
	}
	;

	virtual continuous_post_inters* clone() {
		return new continuous_post_inters(*this);
	}
	;

	/** Return the continuous_set that results from applying time elapse
	 * with the time constraints tcons to the continuous set *cset. */
	virtual continuous::continuous_set_ptr post(const time_constraints& tcons,
			const continuous::continuous_set_const_ptr& cset) const {
		LOGGERSW(DEBUG1,"continuous_post_inters::post","Continuous post with continuous_post_inters");

		typedef continuous::ode_affine_dynamics<scalar_type> ode_aff_dyn;

		continuous::continuous_set_ptr ret_set =
				continuous::continuous_set_ptr();
		if ( continuous::support_function_provider::const_ptr csup =
				boost::dynamic_pointer_cast<const continuous::support_function_provider > (cset)) {
			typename continuous::polyhedron<scalar_type>::const_ptr
					inv_poly_ptr;
			if (tcons.get_invariant()) {
				inv_poly_ptr = boost::dynamic_pointer_cast<
						const continuous::polyhedron<scalar_type> >(
						tcons.get_invariant());
				if (!inv_poly_ptr) {
					std::string tname =
							typeid(*tcons.get_invariant().get()).name();
					throw basic_exception(
							"continuous_post_inters::post cannot handle invariant of type "
									+ tname);
				}
			}

			if (const ode_aff_dyn* dp=dynamic_cast<const ode_aff_dyn*>(tcons.get_dynamics().get())) {
				LOGGER(DEBUG3,"continuous_post_inters::post","computing support_function_time_elapse");

				// call time elapse calculation
				typename continuous::support_function::sfm_cont_set<scalar_type>::ptr
						sfm_ptr =
								continuous::support_function::support_function_time_elapse(
										csup, *dp, inv_poly_ptr,
										get_time_horizon(), get_sampling_time(), my_flowpipe_tolerance);

				LOGGER(DEBUG3,"continuous_post_inters::post","intersecting flowpipe ("+to_string(sfm_ptr->get_size())+" steps, "+to_string(sfm_ptr->get_directions().size())+" directions) with invariant");

				// Intersect the result with the invariant
				// The result is not empty because the root of the sfm is not empty.

				// We do not intersect the invariant with the sfm for precise intersection computation.
//				if (inv_poly_ptr) {
//					 sfm_ptr->intersection_with_poly(*inv_poly_ptr);
//				}

				ret_set = sfm_ptr;
			} else {
				//				std::cerr << "offended by dynamics:"
				//						<< tcons.get_dynamics()->get_predicate() << std::endl;
				std::string tname = typeid(*tcons.get_dynamics().get()).name();
				throw std::runtime_error(
						"continuous_post_inters::post cannot handle dynamics of type "
								+ tname);
			}
		} else {
			std::string tname = typeid(*cset.get()).name();
			throw std::runtime_error(
					"continuous_post_inters::post cannot handle continuous_set of type"
							+ tname);
		}

		return ret_set;
	}
	;

private:
	double my_flowpipe_tolerance;
};

}

#endif /* CONTINUOUS_POST_SFM_H_ */
