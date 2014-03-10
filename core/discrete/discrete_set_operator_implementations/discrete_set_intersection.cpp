/*
 * discrete_set_intersection.cpp
 *
 *  Created on: Sep 8, 2009
 *      Author: frehse
 */

#include "core/discrete/discrete_set_operator_implementations/discrete_set_intersection.h"

#include "core/discrete/discrete_set_implementations.h"

namespace discrete {

discrete_set_ptr intersection_operator<discrete_set_stl_set, discrete_set_stl_set>::implement(
		const discrete_set_stl_set* p1, const discrete_set_stl_set* p2) {
	discrete_set_stl_set::ptr res = discrete_set_stl_set::ptr(new discrete_set_stl_set(*p1));
	res->intersection_assign(*p2);
	return res;
}

discrete_set_ptr intersection_operator<discrete_set_stl_set, singleton_set>::implement(
		const discrete_set_stl_set* p1, const singleton_set* p2) {
	discrete_set_stl_set::ptr res = discrete_set_stl_set::ptr(new discrete_set_stl_set(*p1));
	res->intersection_assign(p2->get_object());
	return res;
}

discrete_set_ptr intersection_operator<singleton_set, discrete_set_stl_set>::implement(
		const singleton_set* p1, const discrete_set_stl_set* p2) {
	return intersection_operator<discrete_set_stl_set, singleton_set>::implement(p2, p1);
}

discrete_set_ptr intersection_operator<singleton_set, singleton_set>::implement(
		const singleton_set* p1, const singleton_set* p2) {
	singleton_set::ptr res = singleton_set::ptr(new singleton_set(*p1));
	res->intersection_assign(*p2);
	return res;
}

}
