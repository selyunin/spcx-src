/*
 * finite_hyperbox_operators.h
 *
 *  Created on: Mar 21, 2013
 *      Author: notroot
 */

#ifndef FINITE_HYPERBOX_OPERATORS_H_
#define FINITE_HYPERBOX_OPERATORS_H_

#include "finite_hyperbox.h"
#include "math/vector_utility.h"

namespace continuous {

/** Intersection */
template<typename scalar_type>
finite_hyperbox<scalar_type> compute_intersection(const finite_hyperbox<scalar_type>& X, const finite_hyperbox<scalar_type>& Y) {
	using namespace math;

	if (X.domain()!=Y.domain()) {
		throw std::runtime_error("compute_intersection: missing implementation for finite_hyperbox with different domains");
	}
	const positional_vdomain& dom = X.domain();

	typedef typename finite_hyperbox<scalar_type>::point_type point_type;
	const point_type& cX = X.get_c();
	const point_type& gX = X.get_g();
	const point_type& cY = Y.get_c();
	const point_type& gY = Y.get_g();
	point_type lX = cX-gX;
	point_type uX = cX+gX;
	point_type lY = cY-gY;
	point_type uY = cY+gY;
	point_type l = max(lX,lY);
	point_type u = min(uX,uY);
	point_type c = (l+u)/scalar_type(2);
	point_type g = u-c;
	finite_hyperbox<scalar_type> Z(c,g,dom);
	return Z;
}

}

#endif /* FINITE_HYPERBOX_OPERATORS_H_ */
