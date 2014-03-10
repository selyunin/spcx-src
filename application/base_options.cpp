/*
 * base_options.cpp
 *
 *  Created on: Dec 30, 2009
 *      Author: frehse
 */

#include "base_options.h"

#include "global/global_types.h"
#include "utility/logger.h"

#include "math/numeric/approx_comparator.h"
#include "math/scalar_types/scalar_with_infinity.h"
#include "math/ode_solving/ode_solver.h"
#include "core/scenarios/reachability_scenario.h"
#include "core/scenarios/scenario_chooser.h"
#include "math/scalar_types/rational.h"

namespace options {

void add_base_options() {
	options::options_processor::config.add_options()("verbosity,v",
			boost::program_options::value<std::string>(),
			"set verbosity level (0,l,m,h,d,D1,D2,D3,D4,D5,D6,D7).");
	options::options_processor::config.add_options()("rel-err",
			boost::program_options::value<std::string>(),
			"Set the relative error for floating point computations (1.0e-12).");
	options::options_processor::config.add_options()("abs-err",
			boost::program_options::value<std::string>(),
			"Set the absolute error for floating point computations (1.0e-22).");
	options::options_processor::config.add_options()("ode-rel-tol",
			boost::program_options::value<std::string>(),
			"Set the relative tolerance for ode solvers (1.0e-6).");
	options::options_processor::config.add_options()("ode-abs-tol",
			boost::program_options::value<std::string>(),
			"Set the absolute tolerance for ode solvers (1.0e-9).");
}

void set_verbosity(options::options_processor::variables_map& vmap) {
	std::string vlevel;
	if (options::options_processor::get_string_option(vmap,"verbosity",vlevel)) {
		if (vlevel.size() == 0)
			throw basic_exception("Unrecognized verbosity level \"" + vlevel + "\"");

		char c = vlevel[0];
		switch (c) {
		case 'l':
			logger::set_active_level(logger_level::LOW);
			break;
		case 'm':
			logger::set_active_level(logger_level::MEDIUM);
			break;
		case 'h':
			logger::set_active_level(logger_level::HIGH);
			break;
		case 'd':
			logger::set_active_level(logger_level::DEBUG);
			break;
		case 'D':
			if (vlevel == "D")
				logger::set_active_level(logger_level::DEBUG);
			if (vlevel == "D1")
				logger::set_active_level(logger_level::DEBUG1);
			else if (vlevel == "D2")
				logger::set_active_level(logger_level::DEBUG2);
			else if (vlevel == "D3")
				logger::set_active_level(logger_level::DEBUG3);
			else if (vlevel == "D4")
				logger::set_active_level(logger_level::DEBUG4);
			else if (vlevel == "D5")
				logger::set_active_level(logger_level::DEBUG5);
			else if (vlevel == "D6")
				logger::set_active_level(logger_level::DEBUG6);
			else if (vlevel == "D7")
				logger::set_active_level(logger_level::DEBUG7);
			break;
		case '0':
			logger::set_active_level(logger_level::OFF);
			break;
		default:
			throw basic_exception("Unrecognized verbosity level \"" + vlevel + "\"");
		}
	}
}

bool check_base_options(options::options_processor::variables_map& vmap) {
	return true;
}

namespace {
template<typename T>
struct rel_error_setter {
	static void implement(std::string c) {
		math::numeric::approx_comparator<T>::set_rel_error_bound(
				from_string<T> (c));
		math::numeric::approx_comparator<scalar_with_infinity<T> >::set_rel_error_bound(
				scalar_with_infinity<T>(from_string<T> (c)));
	}
	;
};
template<typename T>
struct abs_error_setter {
	static void implement(std::string c) {
		math::numeric::approx_comparator<T>::set_abs_error_bound(
				from_string<T> (c));
		math::numeric::approx_comparator<scalar_with_infinity<T> >::set_abs_error_bound(
				scalar_with_infinity<T>(from_string<T> (c)));
	}
	;
};
}

bool apply_base_options(options::options_processor::variables_map& vmap) {
	hybrid_automata::reachability_scenario scen =
			hybrid_automata::scenario_chooser::get_scenario();

	std::string s;
	if (options::options_processor::get_string_option(vmap,"rel-err",s)) {
		global_types::coefficient_type_caller<void, std::string,
				rel_error_setter>::call(s, scen.get_number_type());
	}
	if (options::options_processor::get_string_option(vmap,"abs-err",s)) {
		global_types::coefficient_type_caller<void, std::string,
				abs_error_setter>::call(s, scen.get_number_type());
	}

	if (options::options_processor::get_string_option(vmap,"ode-rel-tol",s)) {
		math::ode::ode_defaults::set_rel_tol(from_string<double> (s));
	}
	if (options::options_processor::get_string_option(vmap,"ode-abs-tol",s)) {
		math::ode::ode_defaults::set_abs_tol(from_string<double> (s));
	}

	return true;
}

}
