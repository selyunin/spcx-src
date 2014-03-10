/*
 * logger_stopwatch.h
 *
 *  Created on: Oct 17, 2010
 *      Author: frehse
 */

#ifndef LOGGER_STOPWATCH_H_
#define LOGGER_STOPWATCH_H_

#include "utility/stopwatch.h"
#include "utility/logger.h"

/** A class for logging the time a process takes.
 *
 * The time is measured from construction to destruction
 * of the object created by the macro.
 *
 * To log the time spent form a certain point in the
 * code until the end of the function use
 * @code
 *   LOGGERSW(XXX,where_msg,what_msg);
 *   ... measured code ...
 * @endcode
 * with XXX being the desired message level,
 * where_msg being a debug level message, and
 * what_msg being at the same time a user level message
 * and a unique identifier for measuring accumulated
 * time spent in this portion of the code.
 *
 * To disable measuring cumulated time, use
 * @code
 *   LOGGERSWNC(XXX,where_msg,what_msg);
 *   ... measured code ...
 * @endcode
 *
 * To measure only the time spent in a certain
 * portion of the code, limit the scope of
 * the stopwatch logger with curly braces:
 * @code
 *   ... some code ...
 *   {
 *      LOGGERSW(XXX,where_msg,what_msg);
 *      ... measured code ...
 *   }
 *   ... the rest of the code ...
 * @endcode
 *
 * To attach a message to a logger_stopwatch,
 * use logger::get_last_id() to obtain the
 * corresponding id:
 * @code
 *   LOGGERSW(XXX,where_msg,what_msg);
 *   logger::logger_id sw_id = logger::get_last_id();
 *   ... measured code ...
 *   LOGGER_ATTACH(XXX,where_msg,"try attaching this",sw_id);
 *   ... more measured code ...
 * @endcode
 * The second message will be attached and
 * finally the timing information will be
 * attached after the second message.
 *
 * @note If an uncaught exception is detected at the time of destruction
 * of the stopwatch, the time is not counted in the statistics
 * and no log output is performed.
 *
 */
class logger_stopwatch {
public:
	typedef enum {
		NO_CUMUL, WITH_CUMUL, ONLY_CUMUL
	} cumulative_mode;
	typedef logger::level level;
	typedef stopwatch::duration_type duration_type;

	/** Construct a stopwatch message with level, where, and what.
	 *
	 * The stopwatch is identified with the what message, and cumulative time is counted
	 * if cumul_mode is WITH_CUMUL or ONLY_CUMUL. For ONLY_CUMUL, only the cumulative
	 * statistics are output.
	 *
	 * @attention l cannot be OFF, since that level is reserved for the active level.
	 */
	logger_stopwatch(level l, const std::string& where_msg,
			const std::string& what_msg,
			cumulative_mode cumul_mode = WITH_CUMUL);

	/** Destructor */
	~logger_stopwatch();

	/** Give an intermediate report on time spent since the last report.
	 *
	 * An optional string s can be attached to the logged message
	 * (in front of the timing info, add leading whitespace yourself). */
	void report_delta(const std::string& s = "");

	/** Report the cumulated time of all logger_stopwatches. */
	static void report_all(logger_level::level l);

private:
	std::string my_where_msg;
	std::string name;
	logger_level::level verbose_level;
	logger::logger_id my_id;
	stopwatch sw;
	cumulative_mode my_cumulative_mode;

	struct data;
	typedef std::map<std::string, data> data_map;
	data_map my_data_map;
	data_map& my_previous_data_map;

	struct data {
		duration_type cumul;
		unsigned int nb_calls;
		data_map local_data;
	};

	static data_map _time_map;
	static int _count;
	static data_map _hierarchy_map;
	static data_map* _current_map;
	static stopwatch _global_watch;

	/** Get the currently active data map. */
	static data_map& current_map();

	/** Set the currently active data map. */
	static void current_map(data_map& m);

	// for sorting the time_map copied to a vector
	typedef std::pair<std::string, data> name_data_pair;
	struct dataCmp {
		bool operator()(const name_data_pair &lhs, const name_data_pair &rhs) {
			return lhs.second.cumul > rhs.second.cumul;
		}
	};

	/** Return the entry corresponding to name in the data map m.
	 *
	 * Creates a new entry if there is no existing one.
	 */
	static data_map::iterator get_or_add_entry(data_map& m,
			const std::string& name);

	/** Update the data map with the duration of the logger_stopwatch identified by its string name
	 *
	 * Returns an iterator to the updated entry. */
	static data_map::const_iterator update_map(data_map& m,
			const std::string& ref, const duration_type& d);

	/** Report the cumulated time of all logger_stopwatches in the data map recursively.
	 *
	 * The output level is used for indentation.
	 * If d is provided, the remaining unaccounted time is reported as well.
	 * */
	static void report(std::ostream& os, const data_map& m, unsigned int output_level = 0, duration_type d = duration_type(), duration_type global_d = duration_type());
};

/** A stopwatch logger that measures cumulated time */
#define LOGGERSW( toto, where, what ) \
    logger_stopwatch logger_stopwatch ## __COUNTER__(logger_level::toto,where,what);

/** A stopwatch logger without cumulating time */
#define LOGGERSWNC( toto, where, what ) \
    logger_stopwatch logger_stopwatch ## __COUNTER__(logger_level::toto,where,what,logger_stopwatch::NO_CUMUL);

/** A stopwatch logger only cumulating time */
#define LOGGERSWOC( toto, where, what ) \
    logger_stopwatch logger_stopwatch ## __COUNTER__(logger_level::toto,where,what,logger_stopwatch::ONLY_CUMUL);

#endif /* LOGGER_STOPWATCH_H_ */
