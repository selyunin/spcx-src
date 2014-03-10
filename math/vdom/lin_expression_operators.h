#ifndef LIN_EXPRESSION_OPERATORS_H_
#define LIN_EXPRESSION_OPERATORS_H_

#include "math/vdom/lin_expression.h"
#include <math.h>


namespace math {

template<typename scalar_type> bool operator==(
		const lin_expression<scalar_type>& v1,
		const lin_expression<scalar_type>& v2) {
	if (!math::numeric::approx_comparator<scalar_type>::is_maybe_equal(
			v1.get_inh_coeff(), v2.get_inh_coeff()))
		return false;
	return v1.get_vdom_vec() == v2.get_vdom_vec();
}
;

template<typename scalar_type> lin_expression<scalar_type> operator+(
		const lin_expression<scalar_type>& l1,
		const lin_expression<scalar_type>& l2) {
	return lin_expression<scalar_type> (l1.get_vdom_vec() + l2.get_vdom_vec(),
			l1.get_inh_coeff() + l2.get_inh_coeff());
}
;

template<typename scalar_type> lin_expression<scalar_type> operator-(
		const lin_expression<scalar_type>& l1,
		const lin_expression<scalar_type>& l2) {
	return lin_expression<scalar_type> (l1.get_vdom_vec() - l2.get_vdom_vec(),
			l1.get_inh_coeff() - l2.get_inh_coeff());
}
;

/** Scalar product. */
template<typename scalar_type> scalar_type scalar_product(const lin_expression<
		scalar_type>& l1, const lin_expression<scalar_type>& l2) {
	return scalar_product(l1.get_vdom_vec(), l2.get_vdom_vec())
			+ (l1.get_inh_coeff() * l2.get_inh_coeff());
}
;

/** Product.
 *
 * Throws if the result is not a linear expression. */
template<typename scalar_type> lin_expression<scalar_type> operator*(
		const lin_expression<scalar_type>& l1,
		const lin_expression<scalar_type>& l2) {
	if (l1.is_homogeneous_coeffs_zero())
		return l1.get_inh_coeff() * l2;
	else if (l2.is_homogeneous_coeffs_zero())
		return l2.get_inh_coeff() * l1;
	else {
		std::stringstream s1;
		s1 << l1;
		std::stringstream s2;
		s2 << l2;
		throw basic_exception("Result of multiplication not a linear expression: (" + s1.str() + ") * ("
				+ s2.str() + ")");
	}
}
;

/** Square root */
/*template<typename scalar_type> lin_expression<scalar_type> sqrt(
		const lin_expression<scalar_type>& l) {

	lin_expression<scalar_type> temp = lin_expression<scalar_type> ((l.get_inh_coeff()));
	return temp;
}
;*/


/** Product of linear expression with scalar */
template<typename scalar_type> lin_expression<scalar_type> operator*(
		const scalar_type& c, const lin_expression<scalar_type>& l) {
	return lin_expression<scalar_type> (c * l.get_vdom_vec(), c
			* l.get_inh_coeff());
}
;

/** Product of scalar with linear expression */
template<typename scalar_type> lin_expression<scalar_type> operator*(
		const lin_expression<scalar_type>& l, const scalar_type& c) {
	return c * l;
}
;

/** Division.
 *
 * Throws if the result is not a linear expression. */
template<typename scalar_type> lin_expression<scalar_type> operator/(
		const lin_expression<scalar_type>& l1,
		const lin_expression<scalar_type>& l2) {
	if (l2.is_homogeneous_coeffs_zero())
		return l1 / l2.get_inh_coeff();
	else {
		std::stringstream s1;
		s1 << l1;
		std::stringstream s2;
		s2 << l2;
		throw basic_exception("Result of division not a linear expression: (" + s1.str() + ") / ("
				+ s2.str() + ")");

	}
}
;

/** Division of scalar by linear expression.
 *
 * Throws if the result is not a linear expression. */
template<typename scalar_type> lin_expression<scalar_type> operator/(
		const scalar_type& l1, const lin_expression<scalar_type>& l2) {
	if (l2.is_homogeneous_coeffs_zero())
		return l1 / l2.get_inh_coeff();
	else {
		std::stringstream s1;
		s1 << l1;
		std::stringstream s2;
		s2 << l2;
		throw basic_exception("Result of division not a linear expression: (" + s1.str() + ") / ("
				+ s2.str() + ")");
	}
}
;

/** Division by scalar */
template<typename scalar_type> lin_expression<scalar_type> operator/(
		const lin_expression<scalar_type>& l, const scalar_type& c) {
	return lin_expression<scalar_type> (l.get_vdom_vec() / c, l.get_inh_coeff()
			/ c);
}
;

/** Negation */
template<typename scalar_type> lin_expression<scalar_type> operator-(
		const lin_expression<scalar_type>& l) {
	return lin_expression<scalar_type> (-l.get_vdom_vec(), -l.get_inh_coeff());
}
;

/** Product.
 *
 * Throws if the result is not a linear expression. */
template<typename scalar_type> lin_expression<scalar_type> pow(
		const lin_expression<scalar_type>& l1,
		const lin_expression<scalar_type>& l2) {
	if (l1.is_homogeneous_coeffs_zero() && l2.is_homogeneous_coeffs_zero()) {
		using std::pow;
		scalar_type a = pow(l1.get_inh_coeff(),l2.get_inh_coeff());
		return lin_expression<scalar_type>(a);
	} else {
		std::stringstream s1;
		s1 << l1;
		std::stringstream s2;
		s2 << l2;
		throw basic_exception("Result of exponentiation not a linear expression: (" + s1.str() + ") * ("
				+ s2.str() + ")");
	}
}
;

}
#endif /*LIN_EXPRESSION_OPERATORS_H_*/
