/*
 * rel_abs_error.h
 *
 *  Created on: Oct 30, 2012
 *      Author: notroot
 */

#ifndef REL_ABS_ERROR_H_
#define REL_ABS_ERROR_H_

#include "comp.h"

namespace math {
namespace numeric {

/** A class for describing absolute and relative error
 *
 * The class is intended to be robust against numerical errors.
 * */
template<typename scalar_type>
class rel_abs_error {
public:
	/** Create error zero */
	rel_abs_error() {
		my_rel = scalar_type(0);
		my_abs = scalar_type(0);
	}
	;

	/** Create an error with relative error r and absolute error a.
	 *
	 * If either is close to or below zero, it snaps to zero. */
	rel_abs_error(const scalar_type& r, const scalar_type& a) {
		// snap r and a to zero if just below
		my_rel = r;
		my_abs = a;
		snap_to_zero();
	}
	;

	/** Returns in tribool whether the scalar interval [x,y] satisfies the error
	 * bounds defined by *this.
	 *
	 * The interval [x,y] is assumed to contain the true value.
	 * The max error is therefore y-x.  */
	tribool is_satisfied(const scalar_type& x, const scalar_type& y) const {
		scalar_type e = abs_error(x, y);
		tribool abs_sat = is_LE(e, my_abs);
		if (definitely(abs_sat))
			return abs_sat;
		else {
			// take the smallest absolute value of x and z
			scalar_type rel_err;
//			if (maybe(is_GE(y, scalar_type(0)) && is_LE(x, scalar_type(0)))) {
//				scalar_type rel_above(0);
//				if (y > my_abs/scalar_type(2))
//					rel_above = rel_error(y - my_abs/scalar_type(2), my_abs/scalar_type(2), y);
//				scalar_type rel_below(0);
//				if (x < my_abs/scalar_type(2))
//					rel_below = rel_error(-my_abs/scalar_type(2) - x, x, -my_abs/scalar_type(2));
//				rel_err = std::max(rel_above, rel_below);
//			} else
				rel_err = rel_error(e, x, y);
			if (rel_err < 0)
				return abs_sat;
			else
				return abs_sat || is_LE(rel_err, my_rel);
		}
	}
	;

	/** Returns satisfaction by an interval itv.
	 *
	 * The interval class must provide the member functions lower() and upper(). */
	template<class interval_type> tribool is_satisfied(const interval_type& itv) const {
		return is_satisfied(itv.lower(), itv.upper());
	}

	/** Create the smallest error bounds to which [x,y] conforms. */
	static rel_abs_error measure_error(const scalar_type& x,
			const scalar_type& y) {
		scalar_type a = abs_error(x, y);
		scalar_type rel_err = rel_error(a, x, y);
		return rel_abs_error(rel_err, a);
	}

	/** Returns satisfaction by an interval itv.
	 *
	 * The interval class must provide the member functions lower() and upper(). */
	template<class interval_type>
	static rel_abs_error measure_error(const interval_type& itv) {
		return measure_error(itv.lower(), itv.upper());
	}

	/** Get the absolute error on a given interval
	 *
	 * Snaps to zero. */
	static scalar_type abs_error(const scalar_type& x, const scalar_type& y) {
		scalar_type e = y - x;
		snap_to_zero(e); // to avoid numerical errors
		return e;
	}

	/** Get the absolute error on a given interval
	 *
	 * Snaps to zero. Returns -1 if x<=0<=y. */
	static scalar_type rel_error(const scalar_type& abs_err,
			const scalar_type& x, const scalar_type& y) {
		if (maybe(is_GE(y, scalar_type(0)) && is_LE(x, scalar_type(0)))) {
			return scalar_type(-1);
		}
		// take the smallest absolute value of x and z
		// not using std::abs so scalar_type doesn't need to implement it
		scalar_type min_abs(x);
		if (min_abs < scalar_type(0))
			min_abs = -min_abs;
		if (y < scalar_type(0)) {
			if (min_abs > -y)
				min_abs = -y;
		} else {
			if (min_abs > y)
				min_abs = y;
		}
		scalar_type b = abs_err / min_abs;
		snap_to_zero(b);
		return b;
	}

	/** Addition */
	rel_abs_error operator+(const rel_abs_error& e) const {
		return rel_abs_error(my_rel + e.my_rel, my_abs + e.my_abs);
	}
	;
	/** Subtraction */
	rel_abs_error operator-(const rel_abs_error& e) const {
		return rel_abs_error(my_rel - e.my_rel, my_abs - e.my_abs);
	}
	;
	/** Scalar multiplication */
	rel_abs_error operator*(const scalar_type& x) const {
		return rel_abs_error(x * my_rel, x * my_abs);
	}
	;
	/** Get the relative error threshold */
	const scalar_type& rel() const {
		return my_rel;
	}
	;
	/** Get the absolute error threshold */
	const scalar_type& abs() const {
		return my_abs;
	}
	;
	/** Less than comparison, returns true if both relative and absolute values are strictly lower */
	tribool operator<(const rel_abs_error& e) const {
		return is_LT(my_rel, e.my_rel) && is_LT(my_abs, e.my_abs);
	}
	;
	/** Less or equal comparison, returns true if both relative and absolute values are lower or equal */
	tribool operator<=(const rel_abs_error& e) const {
		return is_LE(my_rel, e.my_rel) && is_LE(my_abs, e.my_abs);
	}
	;

	template<typename T>
	friend rel_abs_error<T> operator*(const T& x, const rel_abs_error<T>& e);
	template<typename T>
	friend rel_abs_error<T> operator/(const rel_abs_error<T>& e, const T& x);
private:
	/** Replace error by zero if negative or close to zero. */
	void snap_to_zero() {
		snap_to_zero(my_rel);
		snap_to_zero(my_abs);
	}
	/** Replace x by zero if negative or close to zero. */
	static void snap_to_zero(scalar_type& x) {
		if (maybe(is_LT(x, scalar_type(0)))) {
			x == scalar_type(0);
		}
	}
	scalar_type my_rel; // relative error
	scalar_type my_abs; // absolute error
};

/** Scalar multiplication */
template<typename scalar_type>
rel_abs_error<scalar_type> operator*(const scalar_type& x,
		const rel_abs_error<scalar_type>& e) {
	return rel_abs_error<scalar_type> (x * e.my_rel, x * e.my_abs);
}
;
/** Division by a scalar */
template<typename scalar_type>
rel_abs_error<scalar_type> operator/(const rel_abs_error<scalar_type>& e,
		const scalar_type& x) {
	return e * (scalar_type(1) / x);
}
;

/** Output to stream */
template<typename scalar_type>
std::ostream& operator<<(std::ostream& os, const rel_abs_error<scalar_type>& e) {
	os << "e_rel=" << e.rel() << ",e_abs=" << e.abs();
	return os;
}

}
}

#endif /* REL_ABS_ERROR_H_ */
