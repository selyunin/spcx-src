/*
 * context_variable_formatter.h
 *
 *  Created on: Apr 6, 2011
 *      Author: frehse
 */

#ifndef CONTEXT_VARIABLE_FORMATTER_H_
#define CONTEXT_VARIABLE_FORMATTER_H_

#include "variable_formatter.h"
#include "core/predicates/dot_context_lookup.h"

/** A formatter that outputs a name that is unique in the given context. */
class context_variable_formatter: public variable_formatter {
public:
	context_variable_formatter(const std::string& context) :
		my_context(context) {
	}
	;
	virtual ~context_variable_formatter() {
	}
	;
	/** Output the variable to the stream. */
	virtual std::ostream& output(std::ostream& os, const variable& var) {
		os << format(var);
		return os;
	}
	;
	/** Retrieve converted variable name. */
	virtual std::string format(const variable& var) {
		std::string ident = variable::get_name(var.get_id());
		return format(ident);
	}
	;
	/** Converted variable name. */
	virtual std::string format(const std::string& ident) {
		std::string new_str = dot_context::in_context_name(ident, my_context);
		if (new_str.find(".") != std::string::npos) {
			new_str = dot_context::smallest_unique_name(new_str, my_context,
					variable::get_names().begin(), variable::get_names().end());
		}
		return new_str;
	}
	;
	/** Set the context */
	virtual void set_context(const std::string& ident) {
		my_context = ident;
	}
	;
private:
	std::string my_context;
};

#endif /* CONTEXT_VARIABLE_FORMATTER_H_ */
