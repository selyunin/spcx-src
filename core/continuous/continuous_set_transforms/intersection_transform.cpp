#include "core/continuous/continuous_set_transforms/intersection_transform.h"
#include "core/continuous/continuous_set.h"

namespace continuous {

intersection_transform::intersection_transform(const continuous_set_ptr& p) :
	my_set(p) {
}

intersection_transform::~intersection_transform() {
}

const continuous_set_ptr& intersection_transform::get_set() const {
	return my_set;
}

void intersection_transform::get_used_and_modif_variables(variable_id_set& used_vars,
		variable_id_set& modif_vars) const {
	used_vars=my_set->get_variable_ids();
	modif_vars=variable_id_set();
}

void intersection_transform::print(std::ostream& os) const {
	os << my_set;
}


void intersection_transform::accept(
		const_visitor& d) const {
	d.dispatch(this);
}

pre_intersection_transform::pre_intersection_transform(const continuous_set_ptr& p) :
	intersection_transform(p) {
}

pre_intersection_transform::~pre_intersection_transform() {
}

continuous_set_const_ptr pre_intersection_transform::get_relation(continuous_set_const_ptr cset) const {
	return my_set;
}

post_intersection_transform::post_intersection_transform(const continuous_set_ptr& p) :
	intersection_transform(p) {
}

post_intersection_transform::~post_intersection_transform() {
}

continuous_set_const_ptr post_intersection_transform::get_relation(continuous_set_const_ptr cset) const {
	relation_ptr p(my_set->clone());
	p->increase_primedness();
	return my_set;
}

}

