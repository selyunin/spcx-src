/*
 * output_formatter.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef OUTPUT_FORMATTER_H_
#define OUTPUT_FORMATTER_H_

#include <iostream>

#include "math/vdom/variable.h"

/** Forward declarations */
namespace discrete {
class discrete_set;
}
namespace continuous {
class continuous_set;
}
namespace hybrid_automata {
class symbolic_state;
class symbolic_state_collection;
}

namespace io {

/** An output_formatter serializes an object to an output stream in a certain
 * format. */
class output_formatter {
public:
	/** Create an output formatter
	 *
	 * The formatter sends its output to the stream os.
	 */
	explicit output_formatter(std::ostream& os);

	/** Destructor */
	virtual ~output_formatter();

	/** Output a prologue (beginning of file) */
	virtual void prologue();

	/** Output an epilogue (end of file) */
	virtual void epilogue();

	/** Output a discrete set */
	virtual void output(const discrete::discrete_set& d);

	/** Output a continuous set */
	virtual void output(const continuous::continuous_set& c);

	/** Output a symbolic state */
	virtual void output(const hybrid_automata::symbolic_state& sstate);

	/** Output the separator between discrete and continuous set of symbolic state */
	virtual void symbolic_state_element_separator();

	/** Output a symbolic state collection */
	virtual void output(const hybrid_automata::symbolic_state_collection& sstates);

	/** Output the separator between symbolic states of a symbolic state collection */
	virtual void symbolic_state_collection_element_separator();

	/** Get the output stream of the formatter */
	std::ostream& get_os();

	/** Get the output variables of the formatter */
	const variable_id_list& get_output_variables() const;

	/** Set the output variables of the formatter */
	void set_output_variables(variable_id_list vis);

	/** Get the context (component name) of the formatter */
	const std::string& get_context() const;

	/** Set the context (component name) of the formatter */
	void set_context(const std::string& context);

private:
	std::ostream& my_os;
	variable_id_list my_output_variables;
	std::string my_context;
};

}

#endif /* OUTPUT_FORMATTER_H_ */
