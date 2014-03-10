#include "core/continuous/continuous_dynamics/constant_bound_dynamics.h"

#include "core/predicates/valuation_function_tree_utility.h"

namespace continuous {

constant_bound_dynamics::constant_bound_dynamics(const continuous_set_ptr& p) :
	my_set(p) {
}

constant_bound_dynamics::~constant_bound_dynamics() {
}

dynamics_predicate::ptr constant_bound_dynamics::get_predicate() const {
	dynamics_predicate::ptr pred=my_set->get_predicate();
	valuation_functions::increase_primedness(pred);
	return pred;
}

const continuous_set::ptr& constant_bound_dynamics::get_set() const {
	return my_set;
}

variable_id_set constant_bound_dynamics::get_variable_ids() const {
	return my_set->get_variable_ids();
}

variable_id_set constant_bound_dynamics::get_unconstrained_variable_ids() const {
	return variable_id_set();
}

void constant_bound_dynamics::print(std::ostream& os) const {
	// print derivatives as primed, so make a clone with primed vars
	continuous_set_ptr ret_set(my_set->clone());
	ret_set->increase_primedness();
	os << ret_set;
}

void constant_bound_dynamics::accept(const_visitor& d) const {
	d.dispatch(this);
}

}

