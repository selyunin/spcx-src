#ifndef LIN_CONSTRAINT_OPERATORS_H_
#define LIN_CONSTRAINT_OPERATORS_H_

#include "math/vdom/lin_constraint.h"

namespace math {

template<typename scalar_type> bool operator==(const lin_constraint<scalar_type>& v1,
		const lin_constraint<scalar_type>& v2) {
	return v1.get_sign()==v2.get_sign() && v1.get_l()==v2.get_l();
}
;

template<typename scalar_type> lin_constraint<scalar_type> operator+(
		const lin_constraint<scalar_type>& c1, const lin_constraint<scalar_type>& c2) {
	return lin_constraint<scalar_type>(c1.get_l()+c2.get_l(), get_weaker_sign(c1.get_sign(),
			c2.get_sign()));
}
;

template<typename scalar_type> lin_constraint<scalar_type> operator*(const scalar_type& c,
		const lin_constraint<scalar_type>& con) {
	if (c>=scalar_type(0))
		return lin_constraint<scalar_type>(c*con.get_l(), con.get_sign());
	else
		return lin_constraint<scalar_type>(c*con.get_l(), get_flipped_sign(con.get_sign()));
}
;

template<typename scalar_type> lin_constraint<scalar_type> operator*(
		const lin_constraint<scalar_type>& l, const scalar_type& c) {
	return c*l;
}
;

}

#endif /*LIN_CONSTRAINT_OPERATORS_H_*/
