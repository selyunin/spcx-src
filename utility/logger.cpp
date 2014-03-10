/*
 * logger.cpp
 *
 *  Created on: Oct 16, 2010
 *      Author: frehse
 */

#include "logger.h"

#include "utility/basic_exception.h"

std::ostream* logger::log_stream = &std::cout;
logger::level logger::active_log_level = logger_level::MEDIUM;
logger::logger_id logger::last_id = logger::logger_id(2); // this needs to be different from the default attach_id
logger::logger_id logger::last_logged_id = logger::logger_id(1); // this needs to be different from the default attach_id
bool logger::first_message = true;
stopwatch logger::log_sw;

logger::logger(level l, const std::string& where_msg,
		const std::string& what_msg, logger_id attach_id) :
	my_where_msg(where_msg), my_what_msg(what_msg), my_level(l),
			my_id(new_id()) {
	if (l == logger_level::OFF)
		throw basic_exception(
				"logger::logger: log level OFF is reserved for the active level. A message cannot have log level OFF.");
	// if attaching, keep the id
	if (attach_id == last_logged_id)
		my_id = attach_id;
	// copy format of log stream to buffer stream
	if (log_stream)
		my_output.copyfmt(*log_stream);
}

logger::~logger() {
	if (!std::uncaught_exception())
	if (log_stream && my_level <= get_active_level()) {

		std::string msg = my_what_msg + my_output.str();

		if (!msg.empty()) {
			bool attach_to_previous = (my_id == last_logged_id);
			// make a linebreak unless we can attach the message
			// to the message with the requested id
			if (!attach_to_previous && !first_message) {
				*log_stream << std::endl;
			}

			// write the header
			std::string head = header();
			if (!attach_to_previous && !head.empty()) {
				*log_stream << head;
				if (!my_output.str().empty())
					*log_stream << std::endl;
			}

			// swallow one line break of message
			if (!attach_to_previous) {
				char lastchar = *msg.rbegin();
				//				if (lastchar != '\n' && lastchar != '\r') {
				if (lastchar == '\n')
					msg.erase(msg.size() - 1, 1);
				lastchar = *msg.rbegin();
				if (lastchar == '\r')
					msg.erase(msg.size() - 1, 1);
			}

			// write the message
			*log_stream << msg;

			// remember that a message has been written
			first_message = false;
			last_logged_id = my_id;
		} // message not empty
	} // logstream available and level active
}

std::ostream& logger::buffer() {
	return my_output;
}

void logger::set_stream(std::ostream& os) {
	log_stream = &os;
	first_message = true;
}

void logger::terminate() {
	if (!first_message) {
		*log_stream << std::endl;
	}
}

std::ostream* logger::get_stream() {
	return log_stream;
}

logger::level logger::get_active_level() {
	return active_log_level;
}

const logger::logger_id& logger::get_id() const {
	return my_id;
}

void logger::set_active_level(level l) {
	if (l == logger_level::ALWAYS)
		throw basic_exception(
				"logger::set_level: log level ALWAYS is reserved for messages. You cannot set the active log level to ALWAYS.");
	if (l == logger_level::NEVER)
		throw basic_exception(
				"logger::set_level: log level NEVER is reserved for messages. You cannot set the active log level to NEVER.");
	active_log_level = l;
}

std::string logger::time_stamp() const {
	return log_sw.string_value();
}

std::string logger::header() const {
	std::string msg;
	if (active_log_level >= logger_level::DEBUG) {
		msg = time_stamp();
		msg += " ";
	}
	if (active_log_level >= logger_level::DEBUG6) {
		msg += "(";
		msg += my_where_msg + ") ";
	}
	return msg;
}

logger::logger_id logger::new_id() {
	++last_id;
	return last_id;
}

logger::logger_id logger::get_last_id() {
	return last_logged_id;
}

void logger::copyfmt_to(std::ostream& os) {
	if (log_stream)
		os.copyfmt(*log_stream);
}

