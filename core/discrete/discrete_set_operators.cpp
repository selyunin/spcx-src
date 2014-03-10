#include "core/discrete/discrete_set_operators.h"

#include <iostream>
#include <stdexcept>
#include "utility/dispatching/double_dispatch.h"

#include "core/discrete/singleton_set.h"
#include "core/discrete/discrete_set_stl_set.h"

#include "core/discrete/discrete_set_operator_implementations/discrete_set_intersection.h"
#include "core/discrete/discrete_set_operator_implementations/discrete_set_union.h"

namespace discrete {

discrete_set_ptr compute_intersection(
		const discrete_set_const_ptr& p1, const discrete_set_const_ptr& p2) {
	return dispatching::double_dispatch_tc<discrete_set_ptr, intersection_operator,
			discrete_set, discrete_set_typelist>(p1.get(), p2.get());
}

// Operators that possibly take as arguments two discrete_sets of different derived classes
discrete_set_ptr compute_or_assign_intersection(discrete_set_ptr p1,
		const discrete_set_const_ptr& p2) {
	return compute_intersection(p1, p2);
}

discrete_set_ptr compute_union(
		const discrete_set_const_ptr& p1, const discrete_set_const_ptr& p2) {
	//throw std::runtime_error("missing implementation for compute_or_assign_union");
	//return discrete_set_ptr();
	return dispatching::double_dispatch_tc<discrete_set_ptr, union_operator,
			discrete_set, discrete_set_typelist>(p1.get(), p2.get());
}
discrete_set_ptr compute_or_assign_union(discrete_set_ptr p1,
		const discrete_set_const_ptr& p2) {
	return compute_union(p1,p2);
}

discrete_set_ptr compute_difference(
		const discrete_set_const_ptr& p1, const discrete_set_const_ptr& p2) {
	throw std::runtime_error("missing implementation for compute_or_assign_difference");
	return discrete_set_ptr();
//	return dispatching::double_dispatch<discrete_set_ptr, intersection_operator_wrapper,
//			discrete_set, discrete_set_typelist>(p1.get(), p2.get());
}
discrete_set_ptr compute_or_assign_difference(discrete_set_ptr p1,
		const discrete_set_const_ptr& p2) {
	return compute_difference(p1,p2);
}

}

