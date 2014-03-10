#include "core/continuous/continuous_set_transforms/reset_function_transform.h"
#include "core/continuous/continuous_set.h"
#include "utility/shared_ptr_output.h"
#include "core/predicates/node_print_visitor.h"

namespace continuous {

reset_function_transform::reset_function_transform(variable_id var, const function_ptr& f) :
	reset_transform(var), my_function(f) {
}

reset_function_transform::~reset_function_transform() {
}

reset_function_transform::function_ptr reset_function_transform::get_function() {
	return my_function;
}

reset_function_transform::function_const_ptr reset_function_transform::get_function_const() const {
	return my_function;
}

void reset_function_transform::get_used_and_modif_variables(variable_id_set& used_vars,
		variable_id_set& modif_vars) const {
	used_vars=valuation_functions::get_variable_ids(my_function);
	modif_vars=variable_id_set();
	modif_vars.insert(get_variable_id());
}

void reset_function_transform::print(std::ostream& os) const {
	os << variable(get_variable_id()) << " = " << my_function;
}

void reset_function_transform::accept(
		const_visitor& d) const {
	d.dispatch(this);
}

}


