/*
 * predicate_tree_operations.cpp
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#include "io/common_input/predicate_tree_operations.h"

#include "utility/calc_string.h"
#include "utility/stl_helper_functions.h"
#include "utility/share_tree.h"
#include "parse_type_chooser.h"

std::pair<tree::node::ptr, tree::node::ptr> divide_tree(tree::node::ptr p) {
	if (parse_type_chooser::get_bool() == global_types::STD_BOOL)
		return valuation_functions::share_tree<bool>(p);
	else if (parse_type_chooser::get_bool() == global_types::CALC_STR)
		return valuation_functions::share_tree<calc_string>(p);
	else {
		throw std::runtime_error("unsupported bool type in dividing predicate tree");
		return std::make_pair(tree::node::ptr(), tree::node::ptr());
	}
}
