#include "core/continuous/continuous_dynamics/continuous_dynamics_composition.h"
//#include "../continuous_set.h"

/** Forward declarations */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
continuous_set_ptr compute_intersection(
		const continuous_set_const_ptr& p1, const continuous_set_const_ptr& p2);
}

namespace continuous {

bool parallel_compose(continuous_dynamics::ptr& t_ret,
		const continuous_dynamics::const_ptr& t1, const continuous_dynamics::const_ptr& t2) {
	assert(t1);
	assert(t2);
	t_ret=parallel_compose(t1,t2);
	return t_ret;
}

continuous_dynamics::ptr parallel_compose(const continuous_dynamics::const_ptr& t1,
		const continuous_dynamics::const_ptr& t2) {
	assert(t1);
	assert(t2);
	return dispatching::double_dispatch_tc<continuous_dynamics::ptr, continuous_dynamics_composition_operator,
	continuous_dynamics, continuous_dynamics_typelist>(t1.get(), t2.get());
}

continuous_dynamics::ptr continuous_dynamics_composition_operator<constant_bound_dynamics, constant_bound_dynamics>::implement(const constant_bound_dynamics* t1,
		const constant_bound_dynamics* t2) {
	// Intersect continuous sets
	continuous_set::ptr cset=compute_intersection(t1->get_set(),t2->get_set());
	continuous_dynamics::ptr new_dyn=continuous_dynamics::ptr(new constant_bound_dynamics(cset));
	return new_dyn;
}

continuous_dynamics::ptr continuous_dynamics_composition_operator<relation_dynamics, relation_dynamics>::implement(const relation_dynamics* t1,
		const relation_dynamics* t2) {
	// Intersect continuous sets
	continuous_set::ptr cset=compute_intersection(t1->get_relation(),t2->get_relation());
	continuous_dynamics::ptr new_dyn=continuous_dynamics::ptr(new relation_dynamics(cset));
	return new_dyn;
}

}
