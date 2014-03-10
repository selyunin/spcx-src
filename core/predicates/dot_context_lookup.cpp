/*
 * dot_context_lookup.cpp
 *
 *  Created on: Oct 20, 2010
 *      Author: frehse
 */

#include "dot_context_lookup.h"

#include "utility/basic_exception.h"

namespace dot_context {

bool starts_with_context(const std::string& name, const std::string context) {
	size_t cs = context.size();
	// if name and context are equal, it's ok
	if (name==context)
		return true;
	// otherwise the length has to be that of
	// the context + 1 for the period plus at least one character.
	if (name.size() < cs + 2)
		return false;
	return (name.substr(0, cs + 1) == context + ".");
}

std::string in_context_name(const std::string& name, const std::string context) {
	size_t cs = context.size();
	if (starts_with_context(name, context)) {
		if (name.size() > cs + 1)
			return name.substr(cs + 1);
		else
			return "";
	} else
		return name;
}

std::string context_free_name(const std::string& name) {
	size_t cs = name.find_last_of(".");
	if (cs != std::string::npos)
		if (cs + 1 < name.size())
			return name.substr(cs + 1);
		else
			return std::string();
	else
		return name;
}

void throw_if_not_contains_context_free(const variable_id_set& v1,
		const variable_id_set& v2, std::string prologue, std::string epilogue) {
//		std::cout << "testing: ";
//		print_variable_id_set(std::cout,v1);
//		std::cout << " contains ";
//		print_variable_id_set(std::cout,v2);
//		std::cout << "?" << std::endl << std::flush;

	if (!set_contains(v1, v2)) {
		variable_id_set vis(v2);
		set_difference_assign(vis, v1);
		std::stringstream ss("");
		for (variable_id_set::const_iterator it = vis.begin(); it != vis.end(); ++it) {
			if (it != vis.begin())
				ss << ",";
			ss << dot_context::context_free_name(variable::get_name(*it));
		}
		throw basic_exception(
				prologue + ss.str() + epilogue
						+ "\nNote that all symbols must be declared as formal parameters of a component.");
	}
}


}
