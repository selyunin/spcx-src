/*
 * options.h
 *
 *  Created on: Jun 21, 2009
 *      Author: frehse
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <boost/program_options.hpp>

namespace options {

/** Typedef for clean up functions. */
typedef void (*clean_up_function)();

class options_processor {
public:
	typedef boost::program_options::options_description options_description;
	typedef boost::program_options::variables_map variables_map;

	static int read_command_line(int ac, char* av[]);

	/** Interpret a string as if it was the command line. */
	static int read_command_line_from_string(const std::string& cline);

	/** Read config file.
	 *
	 * Throws a basic_exception if a problem occurs.
	 */
	static void read_config_file(const std::string& file_name);

	/** Executes the clean up functions. */
	static int clean_up();

	static variables_map& get_vmap() {
		return vmap;
	}
	;

	static const options_description& get_visible() {
		return visible;
	}
	;

	/** Obtain a string option value.
	 *
	 * Returns true if the option opt_name is in vmap,
	 * and assigns the unquoted value to opt_val.
	 */
	static bool get_string_option(const variables_map& vmap,
			const std::string& opt_name, std::string& opt_val);

	/** Obtain a vector of string option values.
	 *
	 * Returns true if the option opt_name is in vmap,
	 * and assigns the unquoted values to opt_val.
	 */
	static bool get_string_vector_option(const variables_map& vmap,
			const std::string& opt_name, std::vector<std::string>& opt_val);

	static void define_options();
	static void register_clean_up_function(clean_up_function f);

	/** Output the list of options est in vmap. */
	static void print_options(const variables_map& vmap, std::ostream& os);

	static options_description generic;
	static options_description config;
	static options_description hidden;
private:
	static options_description cmdline_options;
	static options_description config_file_options;
	static options_description visible;

	static variables_map vmap;

	static std::vector<clean_up_function> my_clean_funs;
};

}

#endif /* OPTIONS_H_ */
