#include "core/continuous/continuous_set.h"

//#include "../../utility/dispatching/double_dispatch.h"

/** Forward declarations */
namespace continuous {
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
continuous_set_ptr
compute_intersection(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2);
//continuous_set_ptr
//compute_or_assign_intersection(continuous_set_ptr p1,
//		const continuous_set_const_ptr& p2);
}

#include "math/vdom/variable.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"
#include "core/continuous/continuous_set_transforms/sequence_transform.h"
#include "core/continuous/continuous_set_transforms/relation_transform.h"
#include "core/continuous/continuous_set_transforms/intersection_transform.h"

namespace continuous {

continuous_set::ptr continuous_set::get_ptr() {
	continuous_set::ptr p = boost::enable_shared_from_this<continuous_set>::shared_from_this();
	return p;
}

continuous_set::const_ptr continuous_set::get_const_ptr() const {
	continuous_set::const_ptr p = boost::enable_shared_from_this<continuous_set>::shared_from_this();
	return p;
}

math::tribool continuous_set::is_universe() const {
	return math::indeterminate();
}

math::tribool continuous_set::is_disjoint_from(const continuous_set_const_ptr& ps) const {
	assert(ps);
	//continuous_set_ptr sp(clone()); // create a copy of *this
	//sp = compute_or_assign_intersection(sp, ps);
	continuous_set_ptr sp = compute_intersection(get_const_ptr(), ps);
	return sp->is_empty();
}

math::tribool continuous_set::contains(const continuous_set_const_ptr& ps) const {
	assert(ps);
//	if (ps->is_empty())
//		return true;
//	continuous_set_ptr sp = ps->clone(); // create a copy of ps
//	const continuous_set_const_ptr mp=get_const_ptr();
//	sp = compute_or_assign_difference(sp, mp);
//	return sp->is_empty();

	// try to upcast to poly
//	return false;
//	std::cout << "Check "<< typeid(*this).name() << " contains "<< typeid(*ps.get()).name() << std::endl;

	return containment_test(get_const_ptr(), ps);
}

void continuous_set::project_to_variables(const variable_id_set& id_set) {
	variable_id_set vis = get_variable_ids();
	set_difference_assign(vis, id_set);
	existentially_quantify_variables(vis);
}

}

