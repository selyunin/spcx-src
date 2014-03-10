#include "core/continuous/continuous_set_transforms/constant_bound_time_elapse_transform.h"
#include "core/continuous/continuous_set.h"

namespace continuous {

constant_bound_time_elapse_transform::constant_bound_time_elapse_transform(
		const continuous_set::const_ptr& p) :
	my_set(p) {
}

constant_bound_time_elapse_transform::~constant_bound_time_elapse_transform() {
}

relation_const_ptr constant_bound_time_elapse_transform::get_relation(continuous_set_const_ptr cset) const {
	throw std::runtime_error("no get_relation for constant_bound_time_elapse_transform");
	return relation_const_ptr();
}

continuous_set_const_ptr constant_bound_time_elapse_transform::get_set() const {
	return my_set;
}

void constant_bound_time_elapse_transform::get_used_and_modif_variables(variable_id_set& used_vars,
		variable_id_set& modif_vars) const {
	used_vars=my_set->get_variable_ids();
	modif_vars=used_vars;
}

void constant_bound_time_elapse_transform::print(std::ostream& os) const {
	// Construct a primed version of my_set
	continuous_set_ptr c(my_set->clone());
	c->increase_primedness();
	os << c;
}

void constant_bound_time_elapse_transform::accept(const_visitor& d) const {
	d.dispatch(this);
}

}

