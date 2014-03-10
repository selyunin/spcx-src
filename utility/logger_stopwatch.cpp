/*
 * logger_stopwatch.cpp
 *
 *  Created on: Oct 17, 2010
 *      Author: frehse
 */

#include "logger_stopwatch.h"

#include "stl_helper_functions.h" // for to_string
using namespace std;

map<string, logger_stopwatch::data> logger_stopwatch::_time_map;
int logger_stopwatch::_count = 0;
map<string, logger_stopwatch::data> logger_stopwatch::_hierarchy_map;
logger_stopwatch::data_map* logger_stopwatch::_current_map =
		&logger_stopwatch::_hierarchy_map;
stopwatch logger_stopwatch::_global_watch;

logger_stopwatch::logger_stopwatch(level l, const std::string& where_msg,
		const std::string& what_msg, cumulative_mode cumul_mode) :
my_previous_data_map(current_map()) {
	my_where_msg = where_msg;
	name = what_msg;
	verbose_level = l;
	// log the start of the process
	my_cumulative_mode = cumul_mode;
	if (my_cumulative_mode != ONLY_CUMUL) {
		logger start_logger(verbose_level, where_msg, what_msg + "... ");
		my_id = start_logger.get_id();
	} else
		my_id = 0;
	++_count;

	// set my data map as the current data map
	// so that all subsequent loggers pubish their info there
//	IFLOGGER_VAR(verbose_level) {
		if (my_cumulative_mode != NO_CUMUL) {
			logger_stopwatch::data_map::iterator iter = get_or_add_entry(
					current_map(), name);
			current_map(iter->second.local_data);
		}
//	}
}

logger_stopwatch::~logger_stopwatch() {
	--_count;
	// restore the previous data map
	if (my_cumulative_mode != NO_CUMUL) {
		current_map(my_previous_data_map);
	}

	if (!std::uncaught_exception()) {
		// update info
		double cumul_seconds = -1.0;
		if (my_cumulative_mode == WITH_CUMUL || my_cumulative_mode
				== ONLY_CUMUL) {
			// update the global info
			data_map::const_iterator iter = update_map(_time_map, name,
					sw.value());
			cumul_seconds = stopwatch::seconds(iter->second.cumul);

			// update the local info if it's different from the global one
			if (&current_map() != &_time_map)
				update_map(current_map(), name, sw.value());
		}

		IFLOGGER_VAR(verbose_level) {
			std::string intro_msg;
			// if no other message has been logged, make the intro shorter
			if (logger::get_last_id() == my_id) {
				intro_msg = "";
			} else {
				intro_msg = name + " done after ";
			}
			std::string time_msg = to_string(stopwatch::seconds(sw.value()))
					+ "s";
			IFLOGGER(DEBUG) {
				if (my_cumulative_mode == WITH_CUMUL) {
					time_msg += ", cumul " + to_string(cumul_seconds) + "s";
				}
			}
			if (my_cumulative_mode != ONLY_CUMUL) {
				logger(verbose_level, my_where_msg, intro_msg + time_msg, my_id);
			}
		}
	}
}

void logger_stopwatch::report_all(logger_level::level l) {
	IFLOGGER_VAR(l) {
		stringstream s;

		// Sort the map according to duration
		std::vector<name_data_pair> myvec(_time_map.begin(), _time_map.end());
		std::sort(myvec.begin(), myvec.end(), dataCmp());

		// compute the sum
		duration_type map_sum = _global_watch.value();

		for (std::vector<name_data_pair>::const_iterator iter = myvec.begin(); iter
				!= myvec.end(); ++iter) {
			s << setw(9) << fixed << setprecision(3) << stopwatch::seconds(
					iter->second.cumul) << " s, ";
			s << setw(3) << fixed << setprecision(0) << stopwatch::seconds(
					iter->second.cumul) / stopwatch::seconds(map_sum) * 100
					<< "%, ";
			s << setw(6) << iter->second.nb_calls << " x " << iter->first
					<< std::endl;
		}
		LOGGER_OS(ALWAYS,"logger_stopwatch::report_all")
			<< "Cumulative time spent:\n" << s.str();

		stringstream recursive_s;
		report(recursive_s, _hierarchy_map, 0, _global_watch.value(), _global_watch.value());
		LOGGER_OS(ALWAYS,"logger_stopwatch::report_all")
			<< "Hierarchical Profile:\n" << recursive_s.str();
	}
}

void logger_stopwatch::report(std::ostream& os, const data_map& m,
		unsigned int level, duration_type d_total, duration_type global_d) {
	const double low_percent = 0.005;
	const double low_time = 0.0005;

	if (m.begin() != m.end()) {
		// Sort the map according to duration
		std::vector<name_data_pair> myvec(m.begin(), m.end());
		std::sort(myvec.begin(), myvec.end(), dataCmp());

		// If there is no total and global, compute it
		duration_type map_sum;
		for (std::vector<name_data_pair>::const_iterator iter = myvec.begin(); iter
				!= myvec.end(); ++iter) {
			map_sum += iter->second.cumul;
		}
		if (d_total == duration_type())
			d_total = map_sum;
		if (global_d == duration_type())
			global_d = map_sum;

		// create a nice fill string for this level
		std::string fillstr;
		for (unsigned int i = 0; i < level; ++i) {
			fillstr += "   |";
			if (i + 1 >= level)
				fillstr += "--";
		}

		double total_secs = stopwatch::seconds(d_total);
		double global_secs = stopwatch::seconds(global_d);
		duration_type remaining = d_total;
		double remaining_secs = stopwatch::seconds(remaining);
		unsigned int output_count = 0;
		for (std::vector<name_data_pair>::const_iterator iter = myvec.begin(); iter
				!= myvec.end(); ++iter) {
			double iter_secs = stopwatch::seconds(iter->second.cumul);
			if (iter_secs >= low_time
					&& iter_secs >= low_percent * global_secs) {
				os << fillstr;
				os << setw(3) << fixed << setprecision(0)
						<< iter_secs / total_secs * 100 << "% ";
				os << iter->first << " (" << iter->second.nb_calls << "x, ";
				os << fixed << setprecision(0) << iter_secs / global_secs * 100
						<< "% of global)" << std::endl;
				++output_count;
				remaining -= iter->second.cumul;
				remaining_secs = stopwatch::seconds(remaining);
				report(os, iter->second.local_data, level + 1,
						iter->second.cumul, global_d);
			}
		}
		if (output_count > 0 && remaining_secs >= low_time && remaining_secs
				>= low_percent * global_secs) {
			if (d_total == duration_type())
				os << fillstr << setw(3) << fixed << setprecision(3)
						<< stopwatch::seconds(remaining) << " s ";
			else
				os << fillstr << setw(3) << fixed << setprecision(0)
						<< remaining_secs / total_secs * 100 << "% ";
			os << "..." << " (";
			os << fixed << setprecision(0) << remaining_secs / global_secs * 100
					<< "% of global)" << std::endl;;
		}
	}
}

void logger_stopwatch::report_delta(const std::string& s) {
	// output the time since last report and remember the report time
	logger(
			verbose_level,
			my_where_msg,
			name + s + " spent " + to_string(stopwatch::seconds(sw.delta()))
					+ "s since last report");
}

logger_stopwatch::data_map& logger_stopwatch::current_map() {
	return *_current_map;
}

void logger_stopwatch::current_map(data_map& m) {
	_current_map = &m;
}

logger_stopwatch::data_map::iterator logger_stopwatch::get_or_add_entry(
		data_map& m, const std::string& name) {
	map<string, data>::iterator iter = m.find(name);
	if (iter == m.end()) {
		data d = { duration_type(), 0, data_map() };
		m[name] = d;
		iter = m.find(name);
	}
	return iter;
}

logger_stopwatch::data_map::const_iterator logger_stopwatch::update_map(
		data_map& m, const std::string& name, const duration_type& t) {
	map<string, data>::iterator iter = get_or_add_entry(m, name);
	iter->second.cumul += t;
	++iter->second.nb_calls;

	return iter;
}
