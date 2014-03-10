/*
 * logger.h
 *
 *  Created on: Oct 16, 2010
 *      Author: frehse
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#include <stdio.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <set>
#include "stl_helper_functions.h"
#include "logger_level.h"
#include "stopwatch.h"

/** A class for logging messages to a stream
 *
 * A logging message consists of where,what and a level.
 * At the moment of creation, the message is sent to the
 * stream unless the current level is OFF.
 * The where message is only used for at debug levels.
 *
 * The simplest call for non-performance critical calls is
 * @code
 *   LOGGER(XXX,where,what);
 * @endcode
 * where XXX is the desired message level.
 *
 * To avoid the needless construction of unused long messages,
 * use output to a stream buffer:
 * @code
 *   LOGGER_OS(XXX,where) << what;
 * @endcode
 * In that form, the what message is output on a new line
 * rather than attached to the message header.
 * If output to the stream is chained, the message
 * is output after the entire chain is processed:
 * @code
 *   LOGGER_OS(XXX,where) << not << output << until << now;
 * @endcode
 *
 * You may attach one message to another (without any
 * intermediate header or linebreak) by passing its
 * id:
 * @code
 *   logger::logger_id id=LOGGER(XXX,where1,"message1 ");
 *   ... some code ...
 *   LOGGER_ATTACH(XXX,where2,"continuing message1 ",id);
 *   ... some code ...
 *   LOGGER_ATTACH(XXX,where3,"and still continuing message1.",id);
 * @endcode
 * Attached messages have the same id as the message
 * they're attached to, so they may be chained.
 *
 * If the code in between generates itself messages, the
 * attached message is logged as if it had not been attached.
 * To adapt your message depending on whether the
 * attachment has been broken or not, compare
 * to the id of the last message output:
 * @code
 *   logger::logger_id id=LOGGER(XXX,where,what);
 *   ... some code ...
 *   if (logger::get_last_id() == id)
 *      LOGGER_ATTACH(XXX,where,attached_what,id);
 *   else
 *      LOGGER_ATTACH(XXX,where,unattached_what);
 * @endcode
 *
 * To enable certain portions of the code only if the
 * active log level is XXX, use
 * @code
 *   IFLOGGER(XXX) {
 *      ... code ...
 *   }
 * @endcode
 *
 * If the log level is stored in a variable, use
 * @code
 *   IFLOGGER_VAR(a_level_variable) {
 *      ... code ...
 *   }
 * @endcode
 *
 * To stream an object to a string using the logger format,
 * use
 * @code
 * 	 std::string s = logger::formatted_to_string(obj);
 * @endcode
 * This can be a useful alternative to LOGGER_OS for including objects in log messages:
 * @code
 *   LOGGER(XXX,where,"found object " + logger::formatted_to_string(obj));
 * @endcode
 *
 *
 * Starting from debug level, all messages are time stamped
 * with the time since program start.
 *
 * @note As long as an uncaught exception is detected by the logger, all
 * log output is suppressed.
 *
 * @remark Inspired by
 * http://www.drdobbs.com/cpp/201804215;jsessionid=HAFB1XMDT1LCLQE1GHPCKHWATMY32JVN
 * */
class logger {
public:
	typedef logger_level::level level;
	typedef unsigned int logger_id;

	/** Construct an message with level, where, and what.
	 *
	 * If the logger_id of the last logger is provided, the what message is attached
	 * to the last message without starting a new line and without header.
	 *
	 * @attention l cannot be OFF, since that level is reserved for the active level.
	 */
	logger(level l, const std::string& where_msg, const std::string& what_msg="", logger_id id=logger_id(0));

	/** Destructor */
	~logger();

	/** Get a buffered stream to which messages are sent.
	 *
	 * The buffer will be sent to the logstream at
	 * distruction of the logger. */
	std::ostream& buffer();

	/** Get a unique id for a message */
	const logger_id& get_id() const;

	/** Set stream to which messages are sent.
	 *
	 * Setting to 0 disables all output of messages.
	 *
	 * Note: The caller is responsible for
	 * not destroying the stream until after the last logger. */
	static void set_stream(std::ostream& os);

	/** Write the last endl to the stream. */
	static void terminate();

	/** Get the currently active level. */
	static level get_active_level();

	/** Get the currently active level. */
	static void set_active_level(level l);

	/** Get the id of the last message logged. */
	static logger_id get_last_id();

	/** Copy the logger format to the argument stream */
	static void copyfmt_to(std::ostream& os);

	/** Stream an object to a string using the logger format */
	template <class T>
	static std::string formatted_to_string(const T& t);

	/** Get stream to which messages are sent.
	 */
	static std::ostream* get_stream();


private:
	/** Get a string of the current time since
	 * program start. */
	std::string time_stamp() const;

	/** Get the header (time stamp and where). */
	std::string header() const;

	/** Obtain a unique id for a new logger. */
	static logger_id new_id();

	friend class basic_warning;

	static std::ostream* log_stream;
	static level active_log_level;
	static logger_id last_id;
	static logger_id last_logged_id;
	static stopwatch log_sw;
	static bool first_message;
	std::string my_where_msg;
	std::string my_what_msg;
	std::stringstream my_output;
	level my_level;
	logger_id my_id;
};


template<class T>
std::string logger::formatted_to_string(const T& obj) {
	std::stringstream ss;
	copyfmt_to(ss);
	ss << obj;
	return ss.str();
}

//#define LOGGER( toto, where, what ) \
//	if (logger_level::toto > logger::get_active_level()) ; \
//	else logger(logger_level::toto,where,what);
#define LOGGER( toto, where, what ) \
    ( logger_level::toto > logger::get_active_level() ? 0 \
      : logger(logger_level::toto,where,what).get_id() )

#define LOGGER_OS( toto, where ) \
	if (logger_level::toto > logger::get_active_level()) ; \
	else (logger(logger_level::toto,where).buffer())

#define LOGGER_ATTACH( toto, where, what, id ) \
	if (logger_level::toto > logger::get_active_level()) ; \
	else logger(logger_level::toto,where,what,id);

#define IFLOGGER( toto ) \
	if (logger_level::toto <= logger::get_active_level())

#define IFLOGGER_VAR( toto ) \
	if (toto <= logger::get_active_level())

#endif /* LOGGER_H_ */
