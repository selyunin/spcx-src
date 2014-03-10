/*
 * basic_warning.cpp
 *
 *  Created on: Oct 15, 2010
 *      Author: frehse
 */

#include "basic_warning.h"

#include "logger.h"

std::ostream* basic_warning::warning_stream = &std::cout;
basic_warning::warning_type_set basic_warning::active_warnings =
		basic_warning::all_warnings();

std::map<basic_warning::type,unsigned int> basic_warning::warning_count;

basic_warning::basic_warning(const std::string& where_msg,
		const std::string& what_msg, type t) :
	my_where_msg(where_msg), my_what_msg(what_msg), my_type(t) {
	if (warning_stream) {
		std::string msg = what();
		if (msg != "") {
			// @todo This should be handled by passing the logger as
			// a stream (with some appropriate wrapper deriving from basic_ostream)
			// but for now this will have to do.
			if (is_active()) {
				// add warning to counter
				if (warning_count.find(t)==warning_count.end())
					warning_count[t]=1;
				else
					++warning_count[t];

				// output message
				if (warning_stream != logger::get_stream()) {
					*warning_stream << msg << std::endl;
				} else {
					// to avoid interleaving with logged messages, send
					// a logged message
					logger(logger_level::ALWAYS, where_msg,
							"WARNING (" + type_msg(t) + ") " + my_what_msg);
				}
			}
		}
	}
}

basic_warning::~basic_warning() {
}

bool basic_warning::is_active() const {
	return (active_warnings.find(my_type) != active_warnings.end());
}

std::string basic_warning::what() const {
	if (is_active()) {
		return prepare_msg();
	} else {
		return std::string();
	}
}

const basic_warning::warning_type_set& basic_warning::get_active() {
	return active_warnings;
}

void basic_warning::set_active(const warning_type_set& W) {
	active_warnings = W;
}

void basic_warning::activate(type t) {
	active_warnings.insert(t);
}

void basic_warning::deactivate(type t) {
	active_warnings.erase(t);
}

basic_warning::warning_type_set basic_warning::all_warnings() {
	warning_type_set W;
	for (unsigned int i = 0; i < END_OF_WARNING_TYPE_LIST; ++i) {
		W.insert(type(i));
	}
	return W;
}

basic_warning::warning_type_set basic_warning::no_warnings() {
	warning_type_set W;
	return W;
}

void basic_warning::set_stream(std::ostream& os) {
	warning_stream = &os;
}

std::string basic_warning::prepare_msg() const {
	std::string msg = "WARNING " + type_msg(my_type);
	if (!my_where_msg.empty())
		msg += " in " + my_where_msg;
	if (!my_what_msg.empty())
		msg += ": " + my_what_msg;
	return msg;
}

std::string basic_warning::type_msg(type t) {
	switch (t) {
	case AMBIGUOUS_OUTPUT:
		return "ambiguous output";
	case UNUSUAL_OUTPUT:
		return "unusual output";
	case INCOMPLETE_OUTPUT:
		return "incomplete output";
	case UNUSUAL_INPUT:
		return "unusual input";
	default:
		return "";
	}
}

void basic_warning::output_statistics() {
	if (!warning_count.empty()) {
		std::stringstream ss;
		ss << "ATTENTION: Warnings were issued." << std::endl;
		if (warning_stream != logger::get_stream()) {
			*warning_stream << ss;
		} else {
			logger(logger_level::ALWAYS, "Warning statistics",
					ss.str());
		}
	}
	for (std::map<basic_warning::type, unsigned int>::const_iterator it =
			warning_count.begin(); it != warning_count.end(); ++it) {
		std::stringstream ss;
		ss << "Warnings on " << type_msg(it->first) << ": " << it->second << std::endl;
		if (warning_stream != logger::get_stream()) {
			*warning_stream << ss;
		} else {
			logger(logger_level::ALWAYS, "Warning statistics",
					ss.str());
		}
	}
}
