/*
 * parse_type_chooser.cpp
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#include "io/common_input/parse_type_chooser.h"

parse_type_chooser::type parse_type_chooser::get_bool() {
	return my_bool;
}
parse_type_chooser::type parse_type_chooser::get_number() {
	return my_number;
}
void parse_type_chooser::set_bool(type t) {
	my_bool = t;
}
void parse_type_chooser::set_number(type t) {
	my_number = t;
}

parse_type_chooser::type parse_type_chooser::my_bool = global_types::STD_BOOL;
parse_type_chooser::type parse_type_chooser::my_number = global_types::GMP_RATIONAL;
