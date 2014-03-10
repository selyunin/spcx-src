/*
 * bool_node_creation.cpp
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#include "io/common_input/bool_node_creation.h"

#include "utility/calc_string.h"
#include "utility/calc_string_operators.h"
#include "utility/stl_helper_functions.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "io/common_input/parse_type_chooser.h"

namespace predicate_parser {

template<typename T> tree::node::ptr create_const_node_from_string(const std::string& s) {
	return tree::node::ptr(new valuation_functions::const_node<T>(from_string<T> (s)));
}
;

tree::node::ptr create_bool_const_node(std::string s) {
	switch (parse_type_chooser::get_bool()) {
	case global_types::STD_BOOL:
		return create_const_node_from_string<
				global_types::type_selector<global_types::STD_BOOL>::type> (s);
		break;
	case global_types::CALC_STR:
		return create_const_node_from_string<
				global_types::type_selector<global_types::CALC_STR>::type> (s);
		break;
	default:
		throw std::runtime_error("unsupported bool type in predicate_parser");
		return tree::node::ptr();
	}
}

bool is_bool_const_node(const tree::node::ptr& p) {
	switch (parse_type_chooser::get_bool()) {
	case global_types::STD_BOOL:
		if (dynamic_cast<valuation_functions::const_node<bool>*>(p.get()))
					return true;
		break;
	case global_types::CALC_STR:
		if (dynamic_cast<valuation_functions::const_node<calc_string>*>(p.get()))
					return true;
		break;
	default:
		throw std::runtime_error("unsupported bool type in is_bool_const_node");
		break;
	}
	return false;
}

bool get_bool_const_node_value(const tree::node::ptr& p) {
	switch (parse_type_chooser::get_bool()) {
	case global_types::STD_BOOL:
		if (valuation_functions::const_node<bool>* q = dynamic_cast<valuation_functions::const_node<bool>*>(p.get()))
					return q->my_val;
		break;
	case global_types::CALC_STR:
		if (valuation_functions::const_node<calc_string>* q= dynamic_cast<valuation_functions::const_node<calc_string>*>(p.get()))
					return from_string<bool>(q->my_val.get_my_string());
		break;
	default:
		throw std::runtime_error("unsupported bool type in get_bool_const_node_value");
		break;
	}
	return false;
}

}
