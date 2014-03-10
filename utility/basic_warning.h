/*
 * basic_warning.h
 *
 *  Created on: Oct 15, 2010
 *      Author: frehse
 */

#ifndef BASIC_WARNING_H_
#define BASIC_WARNING_H_

#include <iostream>
#include <stdexcept>
#include <string>
#include <set>
#include "stl_helper_functions.h"

/** A basic class for generating warnings.
 *
 * A warning consists of a string message (where,what) and a type.
 * At the moment of creation, the message is sent to the
 * warning stream if the warning type is active.
 *
 * Specific types of warnings can be activated and deactivated.
 * */
class basic_warning {
public:
	typedef enum {
		MISC, AMBIGUOUS_OUTPUT, UNUSUAL_OUTPUT, INCOMPLETE_OUTPUT, UNUSUAL_INPUT, MISSING_IMPLEMENTATION, END_OF_WARNING_TYPE_LIST
	} type;
	typedef std::set<type> warning_type_set;

	/** Construct an warning with a message and type.
	 *
	 * The default type is MISC. */
	basic_warning(const std::string& where_msg, const std::string& what_msg, type t = MISC);

	/** Virtual destructor */
	virtual ~basic_warning();

	/** Returns a string representation of the warning if
	 * the warning is active. */
	virtual std::string what() const;

	/** Returns the currently active warnings. */
	static const warning_type_set& get_active();

	/** Sets the currently active warnings. */
	static void set_active(const warning_type_set& W);

	/** Activate the warning type t. */
	static void activate(type t);

	/** Deactivate the warning type t. */
	static void deactivate(type t);

	/** The set of all warnings. */
	static warning_type_set all_warnings();

	/** No warnings. */
	static warning_type_set no_warnings();

	/** Set stream to which warnings are sent.
	 *
	 * Setting to 0 disables all output of warnings.
	 *
	 * Note: The class does not adopt the pointer.
	 * It's the caller who is responsible for
	 * destroying the stream. */
	static void set_stream(std::ostream& os);

	/** Output statistics on issued warnings
	 *
	 * Output clears count. */
	static void output_statistics();

private:
	/** Returns true if the warning is active. */
	bool is_active() const;

	/** Returns the message of the warning. */
	std::string prepare_msg() const;

	/** Returns the type of warning as a string. */
	static std::string type_msg(type t);

	static std::ostream* warning_stream;
	static warning_type_set active_warnings;

	/** Counters for statistics */
	static std::map<basic_warning::type,unsigned int> warning_count;

	std::string my_where_msg;
	std::string my_what_msg;
	type my_type;
};


#endif /* BASIC_WARNING_H_ */
