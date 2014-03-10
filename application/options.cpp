/*
 * options.cpp
 *
 *  Created on: Jun 21, 2009
 *      Author: frehse
 */

#include "application/options.h"

#include <iostream>
#include <fstream>
#include <string>

#include "utility/basic_exception.h"
#include "utility/split_command_line_args.h"
#include "utility/stl_helper_functions.h"

//#include <vector>
//#include <algorithm>

//namespace po = boost::program_options;

namespace options {

using namespace std;

options_processor::options_description options_processor::generic =
		options_processor::options_description("Generic options");
options_processor::options_description options_processor::config =
		options_processor::options_description("Options");
options_processor::options_description options_processor::hidden =
		options_processor::options_description("Hidden options");
options_processor::options_description options_processor::visible =
		options_processor::options_description("Allowed options");
options_processor::options_description options_processor::cmdline_options =
		options_processor::options_description("Command line options");
options_processor::options_description options_processor::config_file_options =
		options_processor::options_description("Config file options");
//po::options_description options_processor::cmdline_options;
//po::options_description options_processor::config_file_options;

/** Declare static members. */
options_processor::variables_map options_processor::vmap =
		options_processor::variables_map();

void options_processor::define_options() {
	// Declare a group of options that will be
	// allowed only on command line

	generic.add_options()("version", "print version string")("help,h",
			"produce help message")("config,g", boost::program_options::value<
			std::string>(), "Read options from a given config file.");

	config.add_options()(
			"additional-options",
			boost::program_options::value<std::string>(),
			"Pass a string in quotes that is interpreted as further command line options. Does not supersede other options.");

	config.add_options()("show-option-values", boost::program_options::value<
			std::string>()->implicit_value("yes"),
			"If provided, lists the option values passed to the program.");

	cmdline_options.add(generic).add(config).add(hidden);

	config_file_options.add(config).add(hidden);

	visible.add(generic).add(config);

}

int options_processor::read_command_line(int ac, char* av[]) {
	using namespace boost::program_options;
	try {
		// add an empty positional options argument so that non-option arguments are rejected
		boost::program_options::positional_options_description p;
		p.add("non-option-arguments", 0);

		// The following version mangles up the order

		// Note: By default, if a prefix to a long option name is found, then the library "guesses"
		//       the option name. This results in wrong options if one option is a prefix of another.
		//       The following change turns this off.
		//       See: https://svn.boost.org/trac/boost/ticket/860
		int
				parse_style =
						boost::program_options::command_line_style::default_style
								^ boost::program_options::command_line_style::allow_guessing;
		boost::program_options::store(
				boost::program_options::command_line_parser(ac, av).style(parse_style). options(
						cmdline_options).positional(p).run(), vmap);
		//boost::program_options::notify(vmap);
		// so we execute the options one by one
		//		parsed_options p = parse_command_line(ac, av, cmdline_options);
		//		for (unsigned int i = 0; i < p.options.size(); ++i) {
		//			variables_map vm;
		//			parsed_options p2(p.description);
		//			p2.options.push_back(p.options[i]);
		//			store(p2, vm);
		//			notify(vm);
		//
		//			// Don't know how else to parse the options that don't have an argument
		//			process_options(vm);
		//		}
	} catch (exception& e) {
		std::string s = e.what();
		if (s.find("too many positional options") == std::string::npos) {
			throw basic_exception("Failed to parse command line.", e);
		} else {
			basic_exception
					b(
							"Non-option command line arguments are not allowed. Did you forget a - or -- before an option?");
			throw basic_exception("Failed to parse command line.", b);
		}
		return 1;
	}
	return 0;
}

int options_processor::read_command_line_from_string(const std::string& cline) {
	// the following works only under windows
	//std::vector<string> args = boost::program_options::split_winmain(s);

	try {
		if (!cline.empty()) {
			std::string s = cline;
			// remove escapes and quotes
			replace(s, "\\\"", "\"");
			trim_quotes(s);
			std::vector<std::string> arg_vec = split_command_line_args(s);

			// add an empty positional options argument so that non-option arguments are rejected
			boost::program_options::positional_options_description p;
			p.add("non-option-arguments", 0);

			// parse the vector of arguments
			boost::program_options::store(
					boost::program_options::command_line_parser(arg_vec).options(
							cmdline_options).positional(p).run(), vmap);
		}
	} catch (std::exception& e) {
		std::string s = e.what();
		if (s.find("too many positional options") == std::string::npos) {
			throw basic_exception("Failed to parse additional-options.", e);
		} else {
			basic_exception
					b(
							"Non-option command line arguments are not allowed. Did you forget a - or -- before an option?");
			throw basic_exception("Failed to parse additional-options.", b);
		}
		return 1;
	}
	return 0;
}

void options_processor::read_config_file(const string& file_name) {
	try {
		ifstream ifs(file_name.c_str());
		if (ifs.fail())
			throw basic_exception("Unable to open file \"" + file_name + "\".");
		boost::program_options::store(parse_config_file(ifs,
				config_file_options), vmap);
		//boost::program_options::notify(vmap);
	} catch (basic_exception& e) {
		// simply pass it on
		throw basic_exception("Error while reading config file \""
				+ strip_path(file_name) + "\".", e);
	} catch (exception& e) {
		throw basic_exception(
				"Error while parsing config file \"" + file_name
						+ "\".\nTry `--help' for more information on available options.",
				e);
	}
}

int options_processor::clean_up() {
	try {
		for (unsigned int i = 0; i < my_clean_funs.size(); ++i) {
			(*my_clean_funs[i])();
		}
	} catch (exception& e) {
		throw basic_exception("Failed to clean up options_processor", e);
		return 1;
	}
	return 0;
}
void options_processor::register_clean_up_function(clean_up_function f) {
	my_clean_funs.push_back(f);
}

bool options_processor::get_string_option(const variables_map& vmap,
		const std::string& opt_name, std::string& opt_val) {
	if (vmap.count(opt_name)) {
		opt_val = vmap[opt_name].as<std::string> ();
		trim_quotes(opt_val);
		return true;
	} else
		return false;
}

bool options_processor::get_string_vector_option(const variables_map& vmap,
		const std::string& opt_name, std::vector<std::string>& opt_val) {
	if (vmap.count(opt_name)) {
		opt_val = vmap[opt_name].as<std::vector<std::string> > ();
		for (unsigned int i = 0; i < opt_val.size(); ++i) {
			trim_quotes(opt_val[i]);
		}
		return true;
	} else
		return false;
}

void options_processor::print_options(const variables_map& vmap,
		std::ostream& os) {
	for (options::options_processor::variables_map::const_iterator it =
			vmap.begin(); it != vmap.end(); ++it) {
		std::string name = it->first;
		const boost::any & val = it->second.value();
		if (val.type() == typeid(std::string)) {
			std::string strval = boost::any_cast<std::string>(val);
			os << name << " = " << strval << std::endl;
		} else if (val.type() == typeid(std::vector<std::string>)) {
			std::vector<std::string> vecval = boost::any_cast<std::vector<
					std::string> >(val);
			for (unsigned int i = 0; i < vecval.size(); ++i) {
				os << name << "[" + to_string(i) + "] = " << vecval[i]
						<< std::endl;
			}
		} else {
			os << name << " = <unknown type> " << std::endl;
		}
	}
}

std::vector<clean_up_function> options_processor::my_clean_funs;
}
