/*
 * split_command_line_args.h
 *
 *  Created on: Oct 16, 2010
 *      Author: frehse
 */

#ifndef SPLIT_COMMAND_LINE_ARGS_H_
#define SPLIT_COMMAND_LINE_ARGS_H_

#include <string>
#include <vector>

/** Split a string the same way the command line arguments
 * are split on the current operating system.
 *
 * Throws a basic_exception if a problem is encountered.
 */
std::vector<std::string> split_command_line_args(const std::string& s);

#endif /* SPLIT_COMMAND_LINE_ARGS_H_ */
