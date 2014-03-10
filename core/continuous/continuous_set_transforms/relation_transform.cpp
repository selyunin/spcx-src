#include "core/continuous/continuous_set_transforms/relation_transform.h"
#include "core/continuous/continuous_set.h"

namespace continuous {

relation_transform::relation_transform(const relation_ptr& p) :
	relation(p) {
}

relation_transform::~relation_transform() {
}

relation_const_ptr relation_transform::get_relation(continuous_set_const_ptr cset) const {
	return relation;
}

relation::ptr relation_transform::get_relation() {
	return relation::ptr(relation->clone());
}

relation_const_ptr relation_transform::get_relation_const() const {
	return relation;
}


void relation_transform::get_used_and_modif_variables(variable_id_set& used_vars,
		variable_id_set& modif_vars) const {
	used_vars=relation->get_variable_ids();
	modif_vars=used_vars;
}

void relation_transform::print(std::ostream& os) const {
	os << relation;
}

void relation_transform::accept(
		const_visitor& d) const {
	d.dispatch(this);
}

}


