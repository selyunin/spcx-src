/*
 * numeric_comp.h
 *
 *  Created on: 24 oct. 2013
 *      Author: goran
 */

#ifndef NUMERIC_COMP_H_
#define NUMERIC_COMP_H_

#include "tribool.h"

namespace plif {

/** This file defines a how to compare two scalar values numerically. */

/** A class for defining numerical errors and storing default values in static variables
 *
 * Comparisons use the boost tribool class.
 * Error bounds are defined statically for each class.
 *
 * The comparisons are true or false if the same holds for conservative variations
 * within the bounds given by the numerical errors. Otherwise, they are indeterminate.
 *
 * By default, the error is zero, so the exact comparison is used.
 *
 * @author Goran Frehse */
template<typename scalar_type>
class numeric_less {
public:
	/** Returns whether x is strictly less than y */
	tribool operator() (const scalar_type& x, const scalar_type& y) {
		if (augmented(x) < diminished(y)) {
			return true;
		} else if (augmented(y) < diminished(x)) {
			return false;
		} else {
			return indeterminate();
		}
	}
	;

	/** Set the relative error for numerical comparisons */
	static void set_rel_err(const scalar_type& rel_err) {
		my_rel_err_fact_pos = scalar_type(1) + rel_err;
		my_rel_err_fact_neg = scalar_type(1) - rel_err;
	}
	;

	/** Set the absolute error for numerical comparisons */
	static void set_abs_err(const scalar_type& abs_err) {
		my_abs_err = abs_err;
	}
	;

	/** Returns rel_err+1 */
	static const scalar_type& rel_err_fact_pos() {
		return my_rel_err_fact_pos;
	}

	/** Returns rel_err-1 */
	static const scalar_type& rel_err_fact_neg() {
		return my_rel_err_fact_neg;
	}

	/** Returns abs_err */
	static const scalar_type& abs_err() {
		return my_abs_err;
	}

	/** Get the value augmented by the error bounds */
	static scalar_type augmented(const scalar_type& x) {
		if (x >= scalar_type(0)) {
			return std::max(rel_err_fact_pos() * x, x + abs_err());
		} else {
			return std::max(rel_err_fact_neg() * x, x + abs_err());
		}
	}

	/** Get the value augmented by the error bounds */
	static scalar_type diminished(const scalar_type& x) {
		if (x >= scalar_type(0)) {
			return std::max(rel_err_fact_neg() * x, x - abs_err());
		} else {
			return std::max(rel_err_fact_pos() * x, x - abs_err());
		}
	}

private:
	static scalar_type my_rel_err_fact_pos;
	static scalar_type my_rel_err_fact_neg;
	static scalar_type my_abs_err;

};

/** Initialize error bounds to zero by default */
template<typename scalar_type> scalar_type numeric_less<scalar_type>::my_rel_err_fact_pos(
		1);
template<typename scalar_type> scalar_type numeric_less<scalar_type>::my_rel_err_fact_neg(
		1);
template<typename scalar_type> scalar_type numeric_less<scalar_type>::my_abs_err(
		0);

/** Returns as tribool whether x<y when accounting for numerical errors. */
template<typename scalar_type>
tribool is_LT(const scalar_type& x, const scalar_type& y) {
	return numeric_less<scalar_type>()(x, y);
}

/** Returns true if definitely x<y, even when accounting for numerical errors. */
template<typename scalar_type>
bool definitely_is_LT(const scalar_type& x, const scalar_type& y) {
	return definitely(is_LT(x, y));
}

/** Returns true if definitely x>y, even when accounting for numerical errors. */
template<typename scalar_type>
bool definitely_is_GT(const scalar_type& x, const scalar_type& y) {
	return definitely(is_LT(y, x));
}

/** Returns true if possibly x==y when accounting for numerical errors. */
template<typename scalar_type>
bool is_MEQ(const scalar_type& x, const scalar_type& y) {
	return is_indeterminate(numeric_less<scalar_type>()(x, y));
}

}

#endif /* NUMERIC_COMP_H_ */
