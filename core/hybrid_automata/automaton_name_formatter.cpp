/*
 * automaton_name_formatter.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: frehse
 */

#include "automaton_name_formatter.h"

#include "automaton_cache.h"
#include "hybrid_automaton.h"
#include "core/predicates/dot_context_lookup.h"
#include "utility/basic_exception.h"

namespace hybrid_automata {

const int automaton_name_formatter::automaton_name_formatter_index =
		std::ios_base::xalloc();

automaton_name_formatter::automaton_name_formatter(const std::string& name,
		bool override) :
	my_name(name), my_override(override), my_old_formatter(0), my_stream(0) {
}

automaton_name_formatter::automaton_name_formatter(const hybrid_automaton& aut,
		bool override) :
	my_name(aut.get_name()), my_override(override), my_old_formatter(0),
			my_stream(0) {
}

automaton_name_formatter::~automaton_name_formatter() {
	const int index = automaton_name_formatter::automaton_name_formatter_index;

	/** Set the formatter of the stream to the previous value
	 * (likely 0, but that's ok) */
	if (my_stream) {
		my_stream->pword(index) = my_old_formatter;
	}
}

std::ostream& automaton_name_formatter::output(std::ostream& os,
		const std::string& aut) {
	os << aut;
	return os;
}

std::ostream& operator<<(std::ostream& os, const automaton_name_formatter& form) {
	const int index = automaton_name_formatter::automaton_name_formatter_index;

	// don't override if a formatter is already set
	if (os.pword(index) == 0 || form.my_override) {
		if (form.my_stream && form.my_stream != &os)
			throw basic_exception(
					"attempt to apply automaton_name_formatter to two different streams");

		automaton_name_formatter* form_nonconst =
				&const_cast<automaton_name_formatter&> (form);
		/** Remember the old formatter associated with this stream */

		form_nonconst->my_old_formatter
				= reinterpret_cast<automaton_name_formatter*> (os.pword(index));
		form_nonconst->my_stream = &os;

		/** Set form as the new formatter for this stream. */
		os.pword(index) = form_nonconst;
	}

	if (form.my_name != "") {
		/** use_form will be the first formatter passed to the stream.
		 * the automaton name that is output is the one from the last formatter
		 * passed to the stream (the current one).
		 */
		automaton_name_formatter* use_form =
				reinterpret_cast<automaton_name_formatter*> (os.pword(index));
		return use_form->output(os, form.my_name);
	} else
		return os;
}

std::string context_automaton_name_formatter::format(const std::string& name) const {
	//std::cout << "context_automaton_name_formatter formatting " << name << " in context " << my_context;
	std::string new_str = dot_context::in_context_name(name, my_context);
	//	replace_first(new_str, my_context + ".", "");
	if (new_str.find(".") != std::string::npos) {
		std::set<std::string> known =
				hybrid_automaton_cache::get_automaton_names();
		new_str = dot_context::smallest_unique_name(new_str, my_context,
				known.begin(), known.end());
	}
	//std::cout << " -> " << new_str << std::endl;
	return new_str;
}

std::string context_automaton_name_formatter::get() const {
	return format(my_name);
}

std::ostream& context_automaton_name_formatter::output(std::ostream& os,
		const std::string& aut) {
	os << format(aut);
	return os;
}

void context_automaton_name_formatter::set_context(const std::string& ident) {
	my_context = ident;
}

}

