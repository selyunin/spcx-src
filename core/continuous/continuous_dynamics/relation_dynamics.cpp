#include "core/continuous/continuous_dynamics/relation_dynamics.h"

namespace continuous {

relation_dynamics::relation_dynamics(const relation::ptr& p) :
	my_relation(p) {
}

relation_dynamics::~relation_dynamics() {
}

dynamics_predicate::ptr relation_dynamics::get_predicate() const {
	dynamics_predicate::ptr pred=my_relation->get_predicate();
	return pred;
}

variable_id_set relation_dynamics::get_variable_ids() const {
	return get_unprimed_variables(my_relation->get_variable_ids());
}

variable_id_set relation_dynamics::get_unconstrained_variable_ids() const {
	return variable_id_set();
}

void relation_dynamics::print(std::ostream& os) const {
	os << my_relation;
}

relation::const_ptr relation_dynamics::get_relation() const {
	return my_relation;
}

void relation_dynamics::set_relation(relation::ptr p){
	my_relation = p;
}

void relation_dynamics::accept(const_visitor& d) const {
	d.dispatch(this);
}

}
