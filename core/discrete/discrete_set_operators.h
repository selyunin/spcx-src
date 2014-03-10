/*
 * discrete_set_operators.h
 *
 *  Created on: Sep 8, 2009
 *      Author: frehse
 */

#ifndef DISCRETE_SET_OPERATORS_H_
#define DISCRETE_SET_OPERATORS_H_

#include "boost/shared_ptr.hpp"
//#include "discrete_set.h"

/** This file declares the (binary) operators on discrete_set.
 * The actual implementation is given in discrete_set_operator_implementations.h. */

namespace discrete {

class discrete_set;
typedef boost::shared_ptr<discrete_set> discrete_set_ptr;
typedef boost::shared_ptr<const discrete_set> discrete_set_const_ptr;

// Operators that possibly take as arguments two discrete_sets of different derived classes

/** Intersection */
discrete_set_ptr compute_intersection(const discrete_set_const_ptr& p1,
		const discrete_set_const_ptr& p2);
discrete_set_ptr compute_or_assign_intersection(discrete_set_ptr p1,
		const discrete_set_const_ptr& p2);

/** Union */
discrete_set_ptr compute_union(const discrete_set_const_ptr& p1,
		const discrete_set_const_ptr& p2);
discrete_set_ptr compute_or_assign_union(discrete_set_ptr p1,
		const discrete_set_const_ptr& p2);

/** Difference */
discrete_set_ptr compute_difference(const discrete_set_const_ptr& p1,
		const discrete_set_const_ptr& p2);
discrete_set_ptr compute_or_assign_difference(discrete_set_ptr p1,
		const discrete_set_const_ptr& p2);

}

#endif /* DISCRETE_SET_OPERATORS_H_ */
