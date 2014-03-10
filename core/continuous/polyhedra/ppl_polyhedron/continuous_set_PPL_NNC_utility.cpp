/*
 * continuous_set_PPL_NNC_utility.cpp
 *
 *  Created on: Sep 1, 2009
 *      Author: frehse
 */

#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC_utility.h"

#include "core/continuous/polyhedra/ppl_polyhedron/convert_to_ppl.h"
#include "io/common_input/predicate_parser.h"

namespace ppl_polyhedron {

continuous_set_PPL_NNC_ptr parse_PPL_NNC(std::string expr) {
	parse_type_chooser::set_bool(global_types::STD_BOOL);
	parse_type_chooser::set_number(global_types::GMP_RATIONAL);
	tree::node::ptr p = predicate_parser::parse_predicate(expr);
	return convert_to_continuous_set_PPL_NNC(p);
}

}

