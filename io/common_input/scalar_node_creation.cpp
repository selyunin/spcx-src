/*
 * scalar_node_creation.cpp
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#include "io/common_input/scalar_node_creation.h"

#include "math/scalar_types/rational.h"
#include "math/matrix.h"
#include "utility/calc_string.h"
#include "utility/calc_string_operators.h"
#include "utility/stl_helper_functions.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "io/common_input/parse_type_chooser.h"
#include "io/common_input/symbol_table.h"

namespace predicate_parser {

template<typename T> tree::node::ptr create_const_node_from_string(const std::string& s) {
	T value = from_string<T> (s);
	return tree::node::ptr(new valuation_functions::const_node<T>(value));
}
;

tree::node::ptr create_scalar_const_node(const std::string& s) {
	switch (parse_type_chooser::get_number()) {
	case global_types::STD_BOOL:
		return create_const_node_from_string<
				global_types::type_selector<global_types::STD_BOOL>::type> (s);
		break;
	case global_types::STD_INT:
		return create_const_node_from_string<
				global_types::type_selector<global_types::STD_INT>::type> (s);
		break;
	case global_types::STD_DOUBLE:
		return create_const_node_from_string<
				global_types::type_selector<global_types::STD_DOUBLE>::type> (s);
		break;
	case global_types::GMP_RATIONAL:
		return create_const_node_from_string<
				global_types::type_selector<global_types::GMP_RATIONAL>::type> (s);
		break;
	case global_types::CALC_STR:
		return create_const_node_from_string<
				global_types::type_selector<global_types::CALC_STR>::type> (s);
		break;
	default:
		throw std::runtime_error("unsupported scalar type in predicate_parser");
		return tree::node::ptr();
	}
}

tree::node::ptr create_matrix_scalar_const_node(const parser::symbol& s) {
	if(parse_type_chooser::get_number() != global_types::GMP_RATIONAL)
		throw std::runtime_error("Only Rational matrix are implemented.");

	math::matrix<Rational> m(s.dim1, s.dim2);

	for(int j = 0; j < s.dim1; ++j){
		for(int i = 0; i < s.dim2; ++i){
			m(j, i) = from_string<Rational>(boost::any_cast<std::string>(s.my_value(j, i)));
		}
	}
	return tree::node::ptr(new valuation_functions::const_node<math::matrix<Rational> >(m));
}

bool is_scalar_const_node(const tree::node::ptr& p) {
	switch (parse_type_chooser::get_number()) {
	case global_types::STD_INT:
		if (dynamic_cast<valuation_functions::const_node<int>*>(p.get()))
					return true;
		break;
	case global_types::STD_DOUBLE:
		if (dynamic_cast<valuation_functions::const_node<double>*>(p.get()))
					return true;
		break;
	case global_types::GMP_RATIONAL:
		if (dynamic_cast<valuation_functions::const_node<Rational>*>(p.get()))
					return true;
		if (dynamic_cast<valuation_functions::const_node<math::matrix<Rational> >*>(p.get()))
					return true;
		break;
	case global_types::CALC_STR:
		if (dynamic_cast<valuation_functions::const_node<calc_string>*>(p.get()))
					return true;
		break;
	default:
		throw std::runtime_error("unsupported scalar type in predicate_parser");
		break;
	}
	return false;
}

}
