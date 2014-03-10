/*
 * automaton_name_formatter.h
 *
 *  Created on: Apr 7, 2011
 *      Author: frehse
 */

#ifndef AUTOMATON_NAME_FORMATTER_H_
#define AUTOMATON_NAME_FORMATTER_H_

#include <string>
#include <iostream>

/** Forward declaration */
namespace hybrid_automata {
class hybrid_automaton;
}

namespace hybrid_automata {
/** An interface for formatting automaton names for output
 *
 * The default behavior is to use the global
 * automaton_name name.
 *
 * The goal is to have a class that turns a automaton name into an
 * adequate string representation. The class is to be used
 * as a manipulator on stream, i.e., once a automaton_name is
 * inserted in the stream, the inserter must be able to find
 * this class (a reference to the manipulator object) and
 * call it.
 *
 * The use should be as follows:
 * std::cout << some_automaton_name_formatter("") << automaton_name_formatter("x");
 *
 * The first constructor sets the formatter for the stream.
 * The second constructor creates a base class formatter,
 * which calls the first formatter to output "x" as needed.
 *
 * @attention A formatter can only be applied to one stream,
 * otherwise it throws.
 * (This is by convenience, so we don't have to track too many
 * streams).
 * The stream must not be destructed before the formatter!
 *
 * @see
 * http://stackoverflow.com/questions/3114460/need-to-make-context-available-to-c-ostream-insertion-operators
 * http://www.nersc.gov/nusers/resources/PDSF/documentation/pgi/pgC++_lib/stdlibug/str_5412.htm
 *
 * @note To avoid the destruction of the formatter object makes
 * all subsequent inserters to the stream fail (to be checked),
 * the formatter carries a pointer to its stream and the previous formatter.
 * In its destructor, it restores the previous formatter on the stream.
 * Of course this creates a problem if the stream is
 * destroyed before the formatter...
 * This could be fixed by registering a callback on destruction:
 * http://www.cplusplus.com/reference/iostream/ios_base/register_callback/
 * */
class automaton_name_formatter {
public:
	/** Create formatter object for outputting the name.
	 *
	 * Overrides previous settings if override true or unspecified. */
	automaton_name_formatter(const std::string& name, bool override = false);

	/** Create formatter object for outputting the name.
	 *
	 * Overrides previous settings if override true or unspecified. */
	automaton_name_formatter(const hybrid_automaton& name, bool override = false);

	/** Virtual destructor for derived classes.
	 */
	virtual ~automaton_name_formatter();

	/** Output the automaton_name to the stream. */
	virtual std::ostream& output(std::ostream& os, const std::string& aut);

protected:
	std::string my_name;
	bool my_override;
	/** Memorize the old formatter for the stream. */
	automaton_name_formatter* my_old_formatter;
	std::ostream* my_stream;

private:
	/** The index at which the reference to the automaton_name_formatter
	 * is stored in the stream. This is a global value obtained
	 * using ios_base::xalloc().
	 */
	static const int automaton_name_formatter_index;

	friend std::ostream& operator<<(std::ostream& os,
			const automaton_name_formatter& form);

	/** Don't use the default constructor. */
	automaton_name_formatter();
};

/** Set the automaton_name_formatter of an output stream and output the name if set.
 * */
std::ostream
& operator<<(std::ostream& os, const automaton_name_formatter& form);

/** A formatter that produces a context-dependent version of the automaton name
 *
 * This formatter by default overrides previous formatters. */
class context_automaton_name_formatter: public automaton_name_formatter {
public:
	/** Create formatter object for outputting the name in the given context. */
	context_automaton_name_formatter(const std::string& name,
			const std::string& context) :
		automaton_name_formatter(name,true), my_context(context) {
	}
	;
	virtual ~context_automaton_name_formatter() {
	}
	;
	/** Retrieve the formatted name within the context given to the formatter. */
	virtual std::string get() const;

	/** Set the context */
	virtual void set_context(const std::string& ident);

	/** Output the automaton_name to the stream. */
	virtual std::ostream& output(std::ostream& os, const std::string& aut);

private:
	/** Format the name within the context given to the formatter. */
	std::string format(const std::string& name) const;

	std::string my_context;
};

}

#endif /* AUTOMATON_NAME_FORMATTER_H_ */
