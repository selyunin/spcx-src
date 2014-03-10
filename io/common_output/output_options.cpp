/*
 * output_options.cpp
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#include "global/global_types.h" // for __float128 operator<<


#include "output_options.h"

#include <stdexcept>
#include <sstream>
#include <vector>
//#include <boost/tr1/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include "utility/basic_exception.h"
#include "utility/stl_helper_functions.h"
#include "math/vdom/variable.h"
#include "core/hybrid_automata/hybrid_automaton.h"
//#include "core/continuous/polyhedra/polyhedron.h"
//#include "core/symbolic_states/symbolic_state.h"
//#include "core/symbolic_states/symbolic_state_collection.h"
//#include "global/global_types.h"
//#include "core/continuous/sfm_poly/sfm_cont_set.h"

#include "io/JVX_format/JVX_formatter.h"
#include "io/GEN_format/GEN_formatter.h"
#include "io/TXT_format/TXT_formatter.h"
#include "io/INTV_format/INTV_formatter.h"

#include "application/options.h"
#include "io/common_input/input_options.h"
#include "io/SX_format/SX_automaton_formatter.h"

#include "core/analysis_algorithms/monitoring/ha_monitor.h" // for accepting timer as output variable

namespace options {

std::vector<variable_id_list> get_output_variables(
		options::options_processor::variables_map& vmap) {
	std::string sys_name = get_system_name(vmap);
	std::vector<variable_id_list> vis;
	std::vector<std::string> out_vars;
	if (options::options_processor::get_string_vector_option(vmap,"output-variables",out_vars)) {
		vis = std::vector<variable_id_list>(out_vars.size());
		for (unsigned int i = 0; i < out_vars.size(); ++i) {
			std::string text = out_vars[i];
			/* First remove the bounds [...,...] since we don't use them */
			delimited_replace(text, "[", "]", "");

			/* Now separate the csv list */
			boost::tokenizer<boost::escaped_list_separator<char> > tok(text);

			/** If the output variables aren't found, an exception is thrown.
			 * Store the exceptions and throw them at the end of the loop
			 * so that all variables are listed. */
			std::list<std::string> unknown_vars;

			for (boost::tokenizer<boost::escaped_list_separator<char> >::iterator
					it = tok.begin(); it != tok.end(); ++it) {
				std::string token=*it;
				boost::trim(token); // remove whitespace
				try {
					if (token==hybrid_automata::ha_monitor<double>::get_timer_name()) {
						vis[i].push_back(hybrid_automata::ha_monitor<double>::get_timer_id());
					} if (token=="SPACETIME_TIME") {
						vis[i].push_back(variable("SPACETIME_TIME").get_id());
					} else {
						std::string
								var_name =
										valuation_functions::variable_node_creator::context_lookup(
												token, sys_name);
						vis[i].push_back(
								variable::get_or_add_variable_id(var_name, 1));
					}
				} catch ( std::exception& e ) {
					unknown_vars.push_back(token);
				}
			}

			if (!unknown_vars.empty()) {
				std::stringstream ss;
				std::copy(unknown_vars.begin(), unknown_vars.end(),
						infix_ostream_iterator<std::string>(ss, ","));
				basic_exception bex("Could not find variable(s) " + ss.str()
						+ " in component " + sys_name);
				throw basic_exception("Could not set output variables as requested.",bex);
			}
		}
	}
	return vis;
}

std::vector<std::ostream*> define_output_stream(
		options::options_processor::variables_map& vmap) {
	std::vector<std::ostream*> fp;
	std::vector<std::string> out_file_name;
	if (options::options_processor::get_string_vector_option(vmap,"output-file",out_file_name)) {
		fp = std::vector<std::ostream*>(out_file_name.size(), 0);
		for (unsigned int i = 0; i < out_file_name.size(); ++i) {
			fp[i] = new std::ofstream(out_file_name[i].c_str());
		}
	} else
		fp = std::vector<std::ostream*>(1, &std::cout);
	return fp;
}

void destroy_output_stream(std::ostream* fp) {
	if (fp != &std::cout)
		delete fp;
}

void destroy_output_streams(std::vector<std::ostream*> fp) {
	for (unsigned int i = 0; i < fp.size(); ++i) {
		destroy_output_stream(fp[i]);
	}
}

void check_output_format_compatible_with_variables(unsigned int index,
		const variable_id_set& vis,
		options::options_processor::variables_map& vmap) {
	// Get output formats for each stream
	std::vector<std::string> oformat;
	if (options::options_processor::get_string_vector_option(vmap,"output-format",oformat)) {
		// ok, take that value
	} else {
//		oformat = std::vector<std::string>(1, "GEN"); // default output format
		oformat = std::vector<std::string>(1, "INTV"); // default output format
	}

	unsigned int c = vis.size();
	if (index < oformat.size()) {
		std::string f = oformat[index];

		if (f == "GEN") { // format GEN,123
			if (c < 2 || c > 4)
				throw basic_exception(
						"GEN output format only supported for 2 or 3 dimensions, not "
								+ to_string(c) + ".");
		} else if (f == "JVX") {
			if (c < 2 || c > 4)
				throw basic_exception(
						"JVX output format only supported for 2 or 3 dimensions, not "
								+ to_string(c) + ".");
		}
	}
}

io::output_formatter* create_output_formatter(const std::string& f,
		std::ostream& os, options::options_processor::variables_map& vmap) {
	if (f == "TXT") {
		return new io::TXT_formatter(os);
		//	} else if (f == "CON") {
		//		my_output_format = CON;
	} else if (f == "GEN") { // format GEN,123
		double eps = 1.0e-5;
		std::string out_err;
		if (options::options_processor::get_string_option(vmap,"output-error",out_err)) {
			eps = from_string<double>(out_err);
		}
		return new io::GEN_formatter(os, eps);
	} else if (f == "JVX") {
			return new io::JVX_formatter(os);
	} else if (f == "JVX_OFFLINE") {
			return new io::JVX_OFFLINE_formatter(os);
	} else if (f == "INTV") {
			return new io::INTV_formatter(os);
	} else if (f == "TRAJ") {
			return new io::TRAJ_formatter(os);
	} else {
		throw basic_exception("unrecognized output format " + f + ".");
	}
	return 0;
}

std::vector<io::output_formatter*> create_output_formatter(
		options::options_processor::variables_map& vmap) {

	// Define output streams
	std::vector<std::ostream*> fps = define_output_stream(vmap);

	// Get output formats for each stream
	std::vector<std::string> oformat;
	if (options::options_processor::get_string_vector_option(vmap,"output-format",oformat)) {
		// ok, take that value
	}
	else {
		// if none is specified, use GEN
		oformat = std::vector<std::string>(1, "GEN"); // default output format
	}
	if (oformat.size() > fps.size()) {
		throw basic_exception(
				"More output formats defined than output files.");
	}

	// Get output variables for each stream
	std::vector<variable_id_list> vis(fps.size());
	vis = get_output_variables(vmap);

	// create a formatter for each output stream
	std::vector<io::output_formatter*> ofs(fps.size(), 0);
	unsigned int active_format = 0;
	for (unsigned int i = 0; i < fps.size(); ++i) {
		// update the output format if there is another one
		// in the list, otherwise the last one is used
		if (i < oformat.size())
			active_format = i;

		ofs[i] = create_output_formatter(oformat[active_format], *fps[i], vmap);
		if (i < vis.size()) {
			ofs[i]->set_output_variables(vis[i]);
		}
		ofs[i]->set_context(get_system_name(vmap));
	}

	return ofs;
}

void destroy_output_formatter(std::vector<io::output_formatter*> of) {
	for (unsigned int i = 0; i < of.size(); ++i) {
		if (of[i]) {
			std::ostream* fp = &of[i]->get_os();
			delete of[i];
			destroy_output_stream(fp);
		}
	}
}

void add_output_options() {
	options::options_processor::config.add_options()(
			"output-file,o",
			boost::program_options::value<std::vector<std::string> >()->composing(),
			"Output file. Several output files can be specified (one per -o flag). They are paired with a corresponding output format and output variables by the order in which they are given. In case an output file has no corresponding output format, the last specified output format is used.");
	options::options_processor::config.add_options()(
			"output-format,f",
			boost::program_options::value<std::vector<std::string> >()->composing(),
			"Output format:\n- TXT : textual\n- CON : constraints in matrix form\n- GEN : vertices in matrix form\n- JVX : JVX format\n- INTV : [min,max] interval on output variables");
	options::options_processor::config.add_options()(
			"output-variables,a",
			boost::program_options::value<std::vector<std::string> >()->composing(),
			"A comma separated list of variables to output. Continuous sets will be projected to these variables before output.");
	options::options_processor::config.add_options()(
			"output-error",
			boost::program_options::value<std::string>(),
			"Error with which the output is produced, a zero value being exact output (default). Larger values may speed up output. The exact meaning depends on the output format. Active only for GEN format.");
	options::options_processor::config.add_options()(
			"output-system-file",
			boost::program_options::value<std::string>(),
			"Output file to which the (flattened) system automaton is written in SX format. ");
	options::options_processor::config.add_options()(
			"output-format-sfm",
			boost::program_options::value<std::string>(),
			"Specifies how sfms and spacetime_flowpipes are output. Choices are 'code' = internal format (default), 'outer' = outer polyhedral approximation of each element, 'hull' = outer polyhedral approximation of the convex hull of the set.");
}

bool check_output_options(options::options_processor::variables_map& vmap) {
	unsigned int nb_output_files = 0;

	std::vector<std::string> ofiles;
	if (options::options_processor::get_string_vector_option(vmap,"output-file",ofiles)) {
		nb_output_files = ofiles.size();
		for (unsigned int i = 0; i < ofiles.size(); ++i) {
			std::ofstream ofs(ofiles[i].c_str());
			if (!ofs.good()) {
				throw basic_exception("Error opening output file " + ofiles[i]
						+ ".");
			}
			ofs.close();
		}
	}
	// check if output variables are in system
	hybrid_automata::hybrid_automaton::ptr sys = get_system(vmap);
	variable_id_set sys_vars = sys->get_variable_ids();

	// allow the monitor timer
	sys_vars.insert(hybrid_automata::ha_monitor<double>::get_timer_id());
	// allow the spacetime time variable
	sys_vars.insert(variable("SPACETIME_TIME").get_id());

	unsigned int nb_output_variable_lists = 0;
	if (vmap.count("output-variables")) {
		std::vector<variable_id_list> varlists = get_output_variables(vmap);
		nb_output_variable_lists = varlists.size();
		for (unsigned int i = 0; i < nb_output_variable_lists; ++i) {
			variable_id_set vis = list_to_set(varlists[i]);
			throw_if_not_contains(sys_vars, vis,
					"Could not select variable(s) ",
					" as output variables because they don't exist in component "
							+ sys->get_name() + ".");

			// check if output format is compatible
			check_output_format_compatible_with_variables(i, vis, vmap);
		}
	}
	// for the files that don't have matching varlist, check with sys variables
	for (unsigned int i = nb_output_variable_lists; i < nb_output_files; ++i) {
		check_output_format_compatible_with_variables(i, sys_vars, vmap);
	}

	std::string o_sys_file;
	if (options::options_processor::get_string_option(vmap,"output-system-file",o_sys_file)) {
		std::ofstream ofs(o_sys_file.c_str());
		if (!ofs.good()) {
			throw basic_exception("Error opening output system file " + o_sys_file
					+ ".");
		}
		ofs.close();
	}
	return true;
}

bool apply_output_options(options::options_processor::variables_map& vmap) {
	std::string o_sfm_format;
	if (options::options_processor::get_string_option(vmap,"output-format-sfm",o_sfm_format)) {
		io::TXT_formatter::set_sfm_format(o_sfm_format);
	}
	return true;
}

bool apply_post_analysis_output_options(
		options::options_processor::variables_map& vmap) {
	// output system if requested
	std::string o_sys_file;
	if (options::options_processor::get_string_option(vmap,"output-system-file",o_sys_file)) {
		std::ofstream ofs(o_sys_file.c_str());
		if (!ofs.good()) {
			throw std::runtime_error("Error opening output system file "
					+ o_sys_file + ".");
		} else {
			hybrid_automata::hybrid_automaton::ptr sys = get_system(vmap);
			io::SX_automaton_formatter pr(ofs);
			pr.file_prologue();
			sys->accept(pr);
			pr.file_epilogue();
		}
		ofs.close();
	}
	return true;
}

}
