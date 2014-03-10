/*
 * output_options.h
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#include <iostream>
#include <fstream>
#include "boost/shared_ptr.hpp"
#include "math/vdom/variable.h"
#include "application/options.h"
#include "io/common_output/output_formatter.h"

namespace options {

std::vector<variable_id_list> get_output_variables(
		options::options_processor::variables_map& vmap);
std::vector<io::output_formatter*> create_output_formatter(
		options::options_processor::variables_map& vmap);
void destroy_output_formatter(std::vector<io::output_formatter*> of);

void add_output_options();
bool check_output_options(options::options_processor::variables_map& vmap);
bool apply_output_options(options::options_processor::variables_map& vmap);
bool apply_post_analysis_output_options(options::options_processor::variables_map& vmap);

}
