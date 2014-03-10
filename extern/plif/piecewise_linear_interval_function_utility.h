/*
 * piecewise_linear_interval_function_utility.h
 *
 *  Created on: Oct 19, 2012
 *      Author: notroot
 */

#ifndef PIECEWISE_LINEAR_INTERVAL_FUNCTION_UTILITY_H_
#define PIECEWISE_LINEAR_INTERVAL_FUNCTION_UTILITY_H_

#include "piecewise_linear_interval_function.h"

#include <fstream>

namespace plif {

/** This file contains utility functions for pretty output etc. */

/** Stream the vertices of a plf
 *
 * Each breakpoint is output.
 * @todo Show also slopes towards infinity */
inline
void plf_vertex_output(std::ostream& os, const piecewise_linear_function& f) {
	piecewise_linear_function::list_of_points_type::const_iterator iter;
	for (iter = f.get_list().begin(); iter != f.get_list().end(); iter++) {
		if (iter != f.get_list().begin()) // don't output the left of the first
			if (iter->get_y() != iter->get_y_left()) {
				os << iter->get_x() << " " << iter->get_y_left() << std::endl;
				os << std::endl; // empty line forces discontinuity
			}
		os << iter->get_x() << " " << iter->get_y() << std::endl;
		if (iter + 1 != f.get_list().end()) // don't output the right of the last
			if (iter->get_y_right() != iter->get_y()) {
				os << std::endl; // empty line forces discontinuity
				os << iter->get_x() << " " << iter->get_y_right() << std::endl;
			}
	}
}

/** Plot a plf using a system call to graph
 *
 * graph is a standard Linux tool from the plotutils package.
 * The format can be X (window), ps (postscript file), gif (gif file).
 * The filename file_name determines the file name of the output file, whose extension
 * corresponds to the format. It also creates an intermediate file, whose name
 * is file_name concatenated with ".txt".
 *
 * Optionally, additional arguments can be passed to graph by passing them in opt_string.
 * The additional arguments are inserted before the filename in the argument list, and
 * don't need to contain leading or trailing spaces.
 *
 * The function returns the result of the system call to graph.
 * */
inline
int plf_graph(const piecewise_linear_function& f, std::string format,
		std::string file_name, std::string opt_string = "") {
	std::ofstream ofile;

	// write f to file
	std::string fn_f(file_name + ".txt");
	ofile.open(fn_f.c_str());
	plf_vertex_output(ofile, f);
	ofile.close();

	// plot both
	std::string epilog;
	if (format != "X") {
		epilog = "> " + file_name + "." + format;
	}

	std::string command_line = "graph -T" + format + " -C -B " + opt_string
			+ " " + fn_f + " " + epilog;
	int res = std::system(command_line.c_str());
	return res;
}

/** Plot vertices of a plif as a closed polygon
 * */
inline
void plif_vertex_output_closed(std::ostream& os,
		const piecewise_linear_interval_function& h) {

	for (piecewise_linear_function::list_of_points_type::const_iterator iter =
			h.get_upper().get_list().begin();
			iter != h.get_upper().get_list().end(); iter++) {
		if (iter != h.get_upper().get_list().begin()) // don't output the left of the first
			if (iter->get_y() != iter->get_y_left()) {
				os << iter->get_x() << " " << iter->get_y_left() << std::endl;
			}
		os << iter->get_x() << " " << iter->get_y() << std::endl;
		if (iter + 1 != h.get_upper().get_list().end()) // don't output the right of the last
			if (iter->get_y_right() != iter->get_y()) {
				os << iter->get_x() << " " << iter->get_y_right() << std::endl;
			}
	}
	// plot the lower bound right to left so things go in clockwise order
	for (piecewise_linear_function::list_of_points_type::const_reverse_iterator iter =
			h.get_lower().get_list().rbegin();
			iter != h.get_lower().get_list().rend(); iter++) {
		if (iter != h.get_lower().get_list().rbegin()) // don't output the left of the first
			if (iter->get_y() != iter->get_y_right()) {
				os << iter->get_x() << " " << iter->get_y_right() << std::endl;
			}
		os << iter->get_x() << " " << iter->get_y() << std::endl;
		if (iter + 1 != h.get_lower().get_list().rend()) // don't output the right of the last
			if (iter->get_y_left() != iter->get_y()) {
				os << iter->get_x() << " " << iter->get_y_left() << std::endl;
			}
	}
	// close the polygon by repeating the first point
	piecewise_linear_function::list_of_points_type::const_iterator iter =
			h.get_upper().get_list().begin();
	if (iter != h.get_upper().get_list().end()) {
		os << iter->get_x() << " " << iter->get_y() << std::endl;
	}
}

/** Plot a plif using a system call to graph
 *
 * graph is a standard Linux tool from the plotutils package.
 * The format can be X (window), ps (postscript file), gif (gif file).
 * The filename file_name determines the file name of the output file, whose extension
 * corresponds to the format. It also creates two intermediate files, whose names
 * are file_name concatenated with "_f.txt" and "_g.txt", respectively.
 *
 * Optionally, additional arguments can be passed to graph by passing them in opt_string.
 * The additional arguments are inserted before the first filename in the argument list, and
 * don't need to contain leading or trailing spaces.
 *
 *
 * The function returns the result of the system call to graph.
 * */
inline
int plif_graph(const piecewise_linear_interval_function& h, std::string format,
		std::string file_name, std::string opt_string = "") {
	std::ofstream ofile;

	// write lower to file
	std::string fn_f(file_name + "_f.txt");
	ofile.open(fn_f.c_str());
	plf_vertex_output(ofile, h.get_lower());
	ofile.close();
	// write upper to file
	std::string fn_g(file_name + "_g.txt");
	ofile.open(fn_g.c_str());
	plf_vertex_output(ofile, h.get_upper());
	ofile.close();

	// plot both
	std::string epilog;
	if (format != "X") {
		epilog = "> " + file_name + "." + format;
	}
	//  -S 4
	std::string command_line = "graph -T" + format + " -C -B " + opt_string
			+ " -m1 -S 3 " + fn_f + " -s -m2 -S 4 " + fn_g + " " + epilog;
	int res = std::system(command_line.c_str());
	return res;
}

}

#endif /* PIECEWISE_LINEAR_INTERVAL_FUNCTION_UTILITY_H_ */
