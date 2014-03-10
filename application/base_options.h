/*
 * base_options.h
 *
 *  Created on: Dec 30, 2009
 *      Author: frehse
 */

#ifndef BASE_OPTIONS_H_
#define BASE_OPTIONS_H_

#include "application/options.h"

namespace options {

void add_base_options();
void set_verbosity(options::options_processor::variables_map& vmap);
bool check_base_options(options::options_processor::variables_map& vmap);
bool apply_base_options(options::options_processor::variables_map& vmap);

}

#endif /* BASE_OPTIONS_H_ */
