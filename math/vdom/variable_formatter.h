/*
 * variable_formatter.h
 *
 *  Created on: Apr 6, 2011
 *      Author: frehse
 */

#ifndef VARIABLE_FORMATTER_H_
#define VARIABLE_FORMATTER_H_

#include <string>
#include <iostream>

/** Forward declaration */
class variable;

/** An interface for formatting variables for output
 *
 * The default behavior is to use the global
 * variable name.
 *
 * The goal is to have a class that turns a variable into an
 * adequate string representation. The class is to be used
 * as a manipulator on stream, i.e., once a variable is
 * inserted in the stream, the inserter must be able to find
 * this class (a reference to the manipulator object) and
 * call it.
 *
 * The use should be as follows:
 * std::cout << some_variable_formatter("toto") << variable("x");
 *
 * @see
 * http://stackoverflow.com/questions/3114460/need-to-make-context-available-to-c-ostream-insertion-operators
 * http://www.nersc.gov/nusers/resources/PDSF/documentation/pgi/pgC++_lib/stdlibug/str_5412.htm
 *
 * @todo add reference to previous formatter so that it can be reinstated when this formatter
 *       gets destroyed
 * */
class variable_formatter {
public:
	virtual ~variable_formatter() {
	}
	;

	/** Output the variable to the stream. */
	virtual std::ostream& output(std::ostream& os, const variable& var);

private:
	/** The index at which the reference to the variable_formatter
	 * is stored in the stream. This is a global value obtained
	 * using ios_base::xalloc().
	 */
	static const int variable_formatter_index;

	friend std::ostream& operator<<(std::ostream& os,
			const variable_formatter& form);
	friend std::ostream& operator<<(std::ostream& os, const variable& x);
};

/** Set the variable_formatter of an output stream. */
std::ostream& operator<<(std::ostream& os, const variable_formatter& form);

/** Stream the name of a variable. */
std::ostream& operator<<(std::ostream& os, const variable& x);

#endif /* VARIABLE_FORMATTER_H_ */
