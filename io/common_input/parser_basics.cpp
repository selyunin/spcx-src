/*
 * parser_basics.cpp
 *
 *  Created on: Mar 1, 2011
 *      Author: frehse
 */

#include "parser_basics.h"

namespace parser {

bool is_number(const std::string& value) {
	return isdigit(value[0]) || value[0] == '-';
}

}
