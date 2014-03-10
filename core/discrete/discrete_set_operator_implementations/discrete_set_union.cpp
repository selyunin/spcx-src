/*
 * discrete_set_union.cpp
 *
 *  Created on: Mar 22, 2011
 *      Author: mgoyal
 */

#include "core/discrete/discrete_set_operator_implementations/discrete_set_union.h"

#include "core/discrete/discrete_set_implementations.h"

namespace discrete {

discrete_set_ptr union_operator<discrete_set_stl_set, discrete_set_stl_set>::implement(
		const discrete_set_stl_set* p1, const discrete_set_stl_set* p2) {
	discrete_set_stl_set::ptr res = discrete_set_stl_set::ptr(new discrete_set_stl_set(*p1));
	//res->union_assign(*p2);
	return res;
}

discrete_set_ptr union_operator<discrete_set_stl_set, singleton_set>::implement(
		const discrete_set_stl_set* p1, const singleton_set* p2) {
	discrete_set_stl_set::ptr res = discrete_set_stl_set::ptr(new discrete_set_stl_set(*p1));
	res->union_assign(p2->get_object());
	return res;
}

discrete_set_ptr union_operator<singleton_set, discrete_set_stl_set>::implement(
		const singleton_set* p1, const discrete_set_stl_set* p2) {
	return union_operator<discrete_set_stl_set, singleton_set>::implement(p2, p1);
}
discrete_set_ptr union_operator<singleton_set, singleton_set>::implement(
		const singleton_set* p1, const singleton_set* p2) {

	discrete_set_stl_set::ptr res = discrete_set_stl_set::ptr(new discrete_set_stl_set());
	res->union_assign(p1->get_object());
	res->union_assign(p2->get_object());
	//singleton_set::ptr res = singleton_set::ptr(new singleton_set(*p1));
	//res->union_assign(*p2);
	return res;
}

}
