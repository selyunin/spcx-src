/*
 * split_command_line_args.cpp
 *
 *  Created on: Oct 16, 2010
 *      Author: frehse
 */

#include "split_command_line_args.h"

#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <fstream>

#include "stl_helper_functions.h" // for to_string
#include "basic_exception.h"

std::vector<std::string> split_command_line_args(const std::string& s) {
	std::string process_id = to_string(getpid());
	std::string tempf1 = "/tmp/split_command_line_args." + process_id;
	std::string tempf2 = "/tmp/split_command_line_args_output." + process_id;

	std::vector<std::string> arg_vec;
	try {
		int code;

		// write a temporary shell script
		std::ofstream ofs(tempf1.c_str());
		if (ofs.fail())
			throw basic_exception("Could not create temporary file \n" + tempf1
					+ "\n.");
		ofs << "#!/bin/sh\nrm -f "+tempf2+"\nuntil [ \"$*\" = \"\" ]\ndo\necho \"$1\" >> "+tempf2+"\nshift\ndone\n";
		ofs.close();

		// execute the shell script
		std::string command = "/bin/sh " + tempf1 + " " + s; // + " > " + tempf2;
		code = system(command.c_str());
		if (code)
			throw basic_exception("Could not execute command \n" + command
					+ "\n");

		// read the output of the shell script
		std::ifstream ifs(tempf2.c_str());
		if (ifs.fail())
			throw basic_exception("Could not read temporary file \"" + tempf2
					+ "\".");
		std::string line;
		while (getline(ifs, line)) {
			arg_vec.push_back(line);
		}
		ifs.close();

		code = system(("rm -f " + tempf1).c_str());
		if (code)
			throw basic_exception("Could not remove temporary file \n" + tempf1
					+ "\n.");
		code = system(("rm -f " + tempf2).c_str());
		if (code)
			throw basic_exception("Could not remove temporary file \n" + tempf2
					+ "\n.");

	} catch (std::exception& e) {
		throw basic_exception(
				"Could not split the following string into command line arguments:\n"
						+ s + "\n", e);
	}

	return arg_vec;
}
