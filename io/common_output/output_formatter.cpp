/*
 * output_formatter.cpp
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#include "output_formatter.h"

#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"

namespace io {

output_formatter::output_formatter(std::ostream& os) : my_os(os) {
}

output_formatter::~output_formatter() {
}
void output_formatter::prologue() {
}
void output_formatter::epilogue() {
}
void output_formatter::output(const discrete::discrete_set& d) {
}
void output_formatter::output(const continuous::continuous_set& c) {
}
void output_formatter::symbolic_state_element_separator() {
}
void output_formatter::output(const hybrid_automata::symbolic_state& sstate) {
	output(*sstate.get_discrete_set());
	symbolic_state_element_separator();
	output(*sstate.get_continuous_set());
}
void output_formatter::symbolic_state_collection_element_separator() {
}
void output_formatter::output(const hybrid_automata::symbolic_state_collection& sstates) {
	for (hybrid_automata::symbolic_state_collection::const_iterator it =
			sstates.begin(); it != sstates.end(); ++it) {
		if (it != sstates.begin())
			symbolic_state_collection_element_separator();
		output(**it);
	}
}

std::ostream& output_formatter::get_os() {
	return my_os;
}

const variable_id_list& output_formatter::get_output_variables() const {
	return my_output_variables;
}
void output_formatter::set_output_variables(variable_id_list vis) {
	my_output_variables=vis;
}

void output_formatter::set_context(const std::string& context) {
	my_context = context;
}
const std::string& output_formatter::get_context() const {
	return my_context;
}

}
