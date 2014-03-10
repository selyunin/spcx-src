/*
 * breakpoint.h
 *
 *  Created on: 24 oct. 2013
 *      Author: goran
 */

#ifndef BREAKPOINT_H_
#define BREAKPOINT_H_

#include <iostream>

#include "numeric_comp.h"

namespace plif {

/**
 * The class point has four member fields x, y value, left limit at x and right limit at x
 */
class breakpoint {
private:
	precision_type x, y_left, y, y_right; /** A breakpoint is to represent breakpoints, by their x co-ordinate, y value and the left and right limits */
public:
	/**
	 * Constructor
	 */
	breakpoint() {
	}
	/**
	 * Constructor
	 */
	breakpoint(precision_type a, precision_type b, precision_type c,
			precision_type d) {
		x = a;
		y_left = b;
		y = c;
		y_right = d;
	}
	/**
	 * Constructor for a continuous breakpoint
	 *
	 * In a continuous breakpoint, left and right limits are identical to the function value.
	 */
	breakpoint(precision_type a, precision_type b) {
		x = a;
		y_left = b;
		y = b;
		y_right = b;
	}
	/**
	 *Functions to set and get values of the member data fields
	 */
	void set_x(precision_type a) {
		x = a;
	}
	void set_y_left(precision_type a) {
		y_left = a;
	}
	void set_y(precision_type a) {
		y = a;
	}
	void set_y_right(precision_type a) {
		y_right = a;
	}
	const precision_type& get_x() const {
		return x;
	}
	const precision_type& get_y_left() const {
		return y_left;
	}
	const precision_type& get_y() const {
		return y;
	}
	const precision_type& get_y_right() const {
		return y_right;
	}
	precision_type& get_x() {
		return x;
	}
	precision_type& get_y_left() {
		return y_left;
	}
	precision_type& get_y() {
		return y;
	}
	precision_type& get_y_right() {
		return y_right;
	}
	/**
	 * Overloading == operator for the class variable, necessary for this class to be vaid class template type for vector
	 */
	bool operator==(const breakpoint& to_check_equal) const {
		return ((x == to_check_equal.get_x())
				&& (y_left == to_check_equal.get_y_left())
				&& (y == to_check_equal.get_y())
				&& (y_right == to_check_equal.get_y_right()));
	}
	/**
	 * Overloading < operator for the class variable, necessary for this class to be vaid class template type for vector
	 * Might want to redefine this overloading later, as per requirement.
	 */
	bool operator<(const breakpoint& to_check_greater) const {
		return (x < to_check_greater.get_x());
	}
	;

	/** Obtain a breakpoint for sup(f(x)) */
	breakpoint sup_at_x() const {
		return breakpoint(x, sup());
	}
	;

	/** Obtain sup at x */
	const precision_type& sup() const {
		return std::max(y, std::max(y_left, y_right));
	}
	;

	/** Obtain sup at x from the right */
	const precision_type& sup_right() const {
		return std::max(y, y_right);
	}
	;

	/** Obtain sup at x from the left */
	const precision_type& sup_left() const {
		return std::max(y, y_left);
	}
	;

	/** Obtain inf at x */
	const precision_type& inf() const {
		return std::min(y, std::min(y_left, y_right));
	}
	;

	/** Obtain inf at x from the right */
	const precision_type& inf_right() const {
		return std::min(y, y_right);
	}
	;

	/** Obtain inf at x from the left */
	const precision_type& inf_left() const {
		return std::min(y, y_left);
	}
	;

	/** Scalar Multiplication */
	breakpoint& operator*=(const precision_type& c) {
		y_left *= c;
		y *= c;
		y_right *= c;
		return *this;
	}

	/** Scalar Addition */
	breakpoint& operator+=(const precision_type& c) {
		y_left += c;
		y += c;
		y_right += c;
		return *this;
	}

};

/** Stream output for breakpoints
 *
 * @author Goran Frehse */
inline std::ostream& operator<<(std::ostream& os, const breakpoint& p) {
	os << "[" << p.get_x() << ";" << p.get_y_left() << "," << p.get_y() << ","
			<< p.get_y_right() << "]";
	return os;
}
;


/** Mirror a breakpoint along the y axis */
inline
breakpoint breakpoint_x_mirror(const breakpoint& b) {
	return breakpoint(-b.get_x(),b.get_y_left(),b.get_y(),b.get_y_right());
}

/** A struct for determining the computing parameters of a piecewise_linear_function */
struct computation_parameters {
	enum rounding_type { up, down };
	rounding_type rounding;
	/** Default values */
	computation_parameters() {
		rounding = up;
	}
};

/** A function for negating computation parameters
 *
 * rounding up is replaced by rounding down.
 */
inline
computation_parameters negate(const computation_parameters& p) {
	computation_parameters q = p;
	if (p.rounding == computation_parameters::up) {
		q.rounding == computation_parameters::down;
	} else if (p.rounding == computation_parameters::down) {
		q.rounding == computation_parameters::up;
	} else {
		throw std::runtime_error("unknown rounding type");
	}
	return q;
};

/** A function for combining different computation parameters
 *
 * Throws if opposing rounding values are used.
 */
inline
computation_parameters combine(const computation_parameters& p1, const computation_parameters& p2) {
	computation_parameters q = p1;
	if (p1.rounding != p2.rounding) {
		throw std::runtime_error("incompatible rounding modes");
	}
	return q;
};

/** Return a breakpoint obtained by linear interpolation between p and q at position x
 *
 * It is assumed that q is to the right of p.
 * Assuming that x is strictly between p and q, the breakpoint is forcibly continuous
 * and its position determined by the right limit at p and the left limit at q.
 *
 * @author Goran Frehse */
inline
breakpoint interpolate(const breakpoint& p,const breakpoint& q,const precision_type& x, const computation_parameters& params = computation_parameters()) {
	assert(q.get_x()>p.get_x());
	if (x==p.get_x()) {
		return p;
	}
	if (x==q.get_x()) {
		return q;
	}
	if (x < p.get_x() || x > q.get_x()) {
		throw std::runtime_error("interpolating outside of breakpoint range");
	}
	// note that we are forcibly to the right of p.x and to the left of q.x, so what counts is the
	// limits, not the values at the actual breakpoints
	const precision_type& py = p.get_y_right();
	const precision_type& qy = q.get_y_left();
	// because of numerical errors, we could be outside the bounding box of the points. avoid this.
	const precision_type& ymax = std::max(py,qy);
	const precision_type& ymin = std::min(py,qy);

	// if numerical problems arise, round up or down
	if (is_MEQ(p.get_x(),q.get_x())) {
		// it makes no sense computing a slope, so round
		if (params.rounding==computation_parameters::up) {
			return breakpoint(x,ymax);
		} else if (params.rounding==computation_parameters::down) {
			return breakpoint(x,ymin);
		} else {
			throw std::runtime_error("unknown rounding parameter");
		}
	}

	precision_type y = py + (qy-py) / (q.get_x()-p.get_x()) * (x - p.get_x());
	y = std::max(y,ymin);
	y = std::min(y,ymax);
	return breakpoint(x,y);
}

}

#endif /* BREAKPOINT_H_ */
