/*
 * variable_formatter.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: frehse
 */

#include "variable_formatter.h"
#include "variable.h"

const int variable_formatter::variable_formatter_index =
		std::ios_base::xalloc();

std::ostream& variable_formatter::output(std::ostream& os, const variable& var) {
	os << var.get_name();
	return os;
}

std::ostream& operator<<(std::ostream& os, const variable_formatter& form) {
	const int index = variable_formatter::variable_formatter_index;

	//	if (os.pword(index) == 0)
	os.pword(index) = &const_cast<variable_formatter&> (form);

	return os;
}

std::ostream& operator<<(std::ostream& os, const variable& x) {
	const int index = variable_formatter::variable_formatter_index;

	// If no formatter is defined, use the default formatter.
	// Otherwise, use the formatter passed by the stream.
	if (os.pword(index) == 0) {
		return variable_formatter().output(os, x);
	} else {
		variable_formatter* form =
				reinterpret_cast<variable_formatter*> (os.pword(index));
		return form->output(os, x);
	}
}
