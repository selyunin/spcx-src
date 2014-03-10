/*
 * misc_operators.h
 *
 *  Created on: Oct 24, 2012
 *      Author: kateja
 */

#ifndef MISC_OPERATORS_H_
#define MISC_OPERATORS_H_

#include <iostream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "interval.h"
#include "geometry.h"

namespace plif {

/** A sign operator */
template <typename T>
int sgn(const T& val) {
    return (T(0) < val) - (val < T(0));
}

/** Max over 3 arguments */
template <typename T> const T& max(const T& x,const T& y,const T& z) {
	return std::max(std::max(x,y),z);
}
/** Min over 3 arguments */
template <typename T> const T& min(const T& x,const T& y,const T& z) {
	return std::min(std::min(x,y),z);
}
/** Max over 4 arguments */
template <typename T> const T& max(const T& x,const T& y,const T& z,const T& u) {
	return std::max(std::max(x,y),std::max(z,u));
}
/** Max over 6 arguments */
template <typename T> const T& max(const T& x,const T& y,const T& z,const T& u, const T& v, const T& w) {
	return std::max(std::max(std::max(x,y),std::max(u,v)),std::max(z,w));
}
/** Max over 8 arguments */
template <typename T> const T& max(const T& x,const T& y,const T& z,const T& u, const T& v, const T& w, const T& p, const T& q) {
	return std::max(max(x,y,z,u),max(v,w,p,q));
}

/** Max of absolute value over 3 arguments */
template<typename T> T max_abs(const T& x, const T& y, const T& z) {
	return std::max(std::max(std::abs(x), std::abs(y)), std::abs(z));
}
/** Max of absolute value over 4 arguments */
template<typename T> T max_abs(const T& x, const T& y, const T& z, const T& u) {
	return std::max(std::max(std::abs(x), std::abs(y)),
			std::max(std::abs(z), std::abs(u)));
}
/** Max of absolute value over 6 arguments */
template<typename T> T max_abs(const T& x, const T& y, const T& z, const T& u,
		const T& v, const T& w) {
	return std::max(
			std::max(std::max(std::abs(x), std::abs(y)),
					std::max(std::abs(u), std::abs(v))),
			std::max(std::abs(z), std::abs(w)));
}
/** Max of absolute value over 8 arguments */
template<typename T> T max_abs(const T& x, const T& y, const T& z, const T& u,
		const T& v, const T& w, const T& p, const T& q) {
	return std::max(max_abs(x, y, z, u), max_abs(v, w, p, q));
}

/** If the line segment through p(right) and q(left) intersects the line y=c, crosses is set to true and the
 * corresponding x value is returned.
 *
 * This gives the intersection between the values of p and q if q.get_x()>=get_x().
 * Returns false if the intersection point lies to the left of p or to the right of q.
 * */
inline precision_type crosses_threshold(tribool& crosses, const breakpoint& p,
		const breakpoint& q, const precision_type& c) {
	const precision_type& x1 = p.get_x();
	const precision_type& y1 = p.get_y_right();
	const precision_type& x2 = q.get_x();
	const precision_type& y2 = q.get_y_left();

	// check if the two points aren't aligned or the line isn't parallel to y=c
	if (is_MEQ(x1, x2) || is_MEQ(y1, y2)) {
		if ((c >= y1 && c <= y2) || (c >= y2 && c <= y1)) {
			crosses = true;
		} else if (is_MEQ(y1, c) || is_MEQ(y2, c)) {
			crosses = indeterminate();
		} else {
			crosses = false;
		}
		return x1;
	}

	precision_type x_c = (c - y1) / (y2 - y1) * (x2 - x1) + x1;

	crosses = (!is_LT(x_c,x1) && !is_LT(x2,x_c));

	// limit x_c to [x1,x2]
	x_c = std::max(x1,x_c);
	x_c = std::min(x2,x_c);
	return x_c;
}
;

}

#endif /* MISC_OPERATORS_H_ */
