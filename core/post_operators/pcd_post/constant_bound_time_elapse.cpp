#include "core/continuous/continuous_set_operator_implementations/compute_transformation.h"
#include "core/post_operators/pcd_post/constant_bound_time_elapse.h"
#include "core/continuous/continuous_dynamics/constant_bound_dynamics.h"
#include "core/hybrid_automata/location.h"
#include "utility/shared_ptr_output.h"
#include "core/predicates/node_print_visitor.h"

/** Forward declarations. */
namespace continuous {
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
continuous_set_ptr compute_intersection(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2);
continuous_set_ptr compute_or_assign_intersection(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2);
}
#include "core/continuous/continuous_set_transforms/constant_bound_time_elapse_transform.h"

namespace hybrid_automata {

unsigned int constant_bound_time_elapse_post::max_refine_count = 1;

continuous::continuous_set_ptr constant_bound_time_elapse_post::post(
		const time_constraints& tcons,
		const continuous::continuous_set_const_ptr& cset) const {

	using namespace continuous;

	continuous_set_const_ptr bounds;
	continuous_set_ptr ret_set = continuous_set_ptr();

	// If the dynamics are already constant bounds, we just need to grab them.
	if (const constant_bound_dynamics* dp=dynamic_cast<const constant_bound_dynamics*>(tcons.get_dynamics().get())) {
		bounds = dp->get_set();

		if (!bounds->is_empty()) {
			// apply the transform
			ret_set = compute_transformation(*cset,
					constant_bound_time_elapse_transform(bounds));

			// Intersect the result with the invariant
			ret_set = compute_or_assign_intersection(ret_set,
					tcons.get_invariant());
		} else {
			ret_set = continuous::continuous_set_ptr(cset->clone());
		}
	} else if (const relation_dynamics* dp=dynamic_cast<const relation_dynamics*>(tcons.get_dynamics().get())) {
		// Iteraterate projection and time elapse until fixpoint is found
		continuous_set_ptr new_set = continuous_set_ptr(
				tcons.get_invariant()->clone());

		// compute if it's the first time or if progress has been made
		unsigned int count = 0;
		while (count <= max_refine_count && (!ret_set || !math::definitely(new_set->contains(ret_set)))) {
			++count;
			ret_set = new_set;

			continuous_set_ptr flow_relation;
			// affine or nonlinear dynamics, so we need to project the flow relation
			flow_relation = compute_intersection(dp->get_relation(), new_set); // intersect with the invariant (to gain accuracy in the projection)

			variable_id_set vis = flow_relation->get_primed_variables(0); // get the unprimed variables = state variables
			flow_relation->existentially_quantify_variables(vis); // remove them, leaving us with the derivatives as primed variables
			flow_relation->decrease_primedness(); // we need the derivatives as unprimed variables
			bounds = flow_relation;

			if (!bounds->is_empty()) {
				// apply the transform
				new_set = compute_transformation(*cset,
						constant_bound_time_elapse_transform(bounds));

				// Intersect the result with the invariant
				new_set = compute_or_assign_intersection(new_set,
						tcons.get_invariant());
			} else {
				new_set = continuous::continuous_set_ptr(cset->clone());
			}
		}
		ret_set = new_set;
	} else {
		std::cerr << "offended by dynamics:"
				<< tcons.get_dynamics()->get_predicate() << std::endl;
		throw std::runtime_error(
				"cannot apply constant_bound_time_elapse_post to these dynamics");
	}

	return ret_set;
}

}
