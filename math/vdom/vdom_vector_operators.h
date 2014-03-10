
#ifndef vdom_VECTOR_OPERATORS_H_
#define vdom_VECTOR_OPERATORS_H_

#include "math/numeric/approx_comparator.h"
#include "math/vector_operators.h"
#include "math/vdom/vdom_vector.h"

namespace math {

/** Compare if v1 and v2 are equal. Uses 
 * math::numeric::approx_comparator<scalar_type>::is_maybe_equal
 * for approximative computations.
 * 
 * @todo: right now the complexity is 2*n log n. Optimize.
 */
template<typename scalar_type> bool operator==(
		const vdom_vector<scalar_type>& v1,
		const vdom_vector<scalar_type>& v2) {
	if (v1.get_index_to_variable_id_map() == v2.get_index_to_variable_id_map()) {
		return v1.get_vector() == v2.get_vector();
	} else {
		/* Check that all the nonzero coefficients in v1 exist in v2 */
		for (typename vdom_vector<scalar_type>::size_type i = 0; i < v1.size(); ++i) {
			if (!math::numeric::approx_comparator<scalar_type>::is_maybe_zero(
					v1[i])) {
				variable_id id = v1.get_id(i);
				if (!v2.has_id(id) || !math::numeric::approx_comparator<
						scalar_type>::is_maybe_equal(v2[v2.get_index(id)],
						v1[i]))
					return false;
			} else {
				variable_id id = v1.get_id(i);
				if (v2.has_id(id) && !math::numeric::approx_comparator<
						scalar_type>::is_maybe_zero(v2[v2.get_index(id)]))
					return false;
			}
		}
		/* Check that all the nonzero coefficients in v2 exist in v1 */
		for (typename vdom_vector<scalar_type>::size_type i = 0; i < v2.size(); ++i) {
			if (!math::numeric::approx_comparator<scalar_type>::is_maybe_zero(
					v2[i])) {
				variable_id id = v2.get_id(i);
				if (!v1.has_id(id) || !math::numeric::approx_comparator<
						scalar_type>::is_maybe_equal(v1[v1.get_index(id)],
						v2[i]))
					return false;
			} else {
				variable_id id = v2.get_id(i);
				if (v1.has_id(id) && !math::numeric::approx_comparator<
						scalar_type>::is_maybe_zero(v1[v1.get_index(id)]))
					return false;
			}
		}
	}
	return true;
}
;

template<typename scalar_type> vdom_vector<scalar_type> operator+(
		const vdom_vector<scalar_type>& v1,
		const vdom_vector<scalar_type>& v2) {
	if (v1.get_index_to_variable_id_map() == v2.get_index_to_variable_id_map()) {
		return vdom_vector<scalar_type> (v1.get_vector() + v2.get_vector(),
				v1.get_index_to_variable_id_map());
	} else {
		index_to_variable_id_map_ptr new_iimap;
		index_type new_dim;
		index_to_index_bimap remap2;
		get_common_map(v1.get_index_to_variable_id_map(),
				v2.get_index_to_variable_id_map(), new_iimap, new_dim, remap2);
		typename vdom_vector<scalar_type>::vector_type new_v(new_dim);
		vdom_vector<scalar_type> new_l(new_iimap);
		for (unsigned int i = 0; i < v1.size(); ++i) {
			new_l[i] = v1[i];
		}
		for (unsigned int i = 0; i < v2.size(); ++i) {
			new_l[remap2.get_map(i)] = new_l[remap2.get_map(i)] + v2[i];
		}
		return new_l;
	}
}
;

template<typename scalar_type> vdom_vector<scalar_type> operator-(
		const vdom_vector<scalar_type>& v1,
		const vdom_vector<scalar_type>& v2) {
	if (v1.get_index_to_variable_id_map() == v2.get_index_to_variable_id_map()) {
		return vdom_vector<scalar_type> (v1.get_vector() - v2.get_vector(),
				v1.get_index_to_variable_id_map());
	} else {
		index_to_variable_id_map_ptr new_iimap;
		index_type new_dim;
		index_to_index_bimap remap2;
		get_common_map(v1.get_index_to_variable_id_map(),
				v2.get_index_to_variable_id_map(), new_iimap, new_dim, remap2);
		typename vdom_vector<scalar_type>::vector_type new_v(new_dim);
		vdom_vector<scalar_type> new_l(new_iimap);
		for (unsigned int i = 0; i < v1.size(); ++i) {
			new_l[i] = v1[i];
		}
		for (unsigned int i = 0; i < v2.size(); ++i) {
			new_l[remap2.get_map(i)] = new_l[remap2.get_map(i)] - v2[i];
		}
		return new_l;
	}
}
;

/** Scalar product.
 *
 * For two vectors \f$ a,b \f$ returns \f$ \sum_{x \in V} a_x b_x \f$, where
 * \f$ V = V_a \cup V_b \f$ and \f$ a_x=0 \f$ for \f$ x \notin V_a \f$,
 * \f$ b_x=0 \f$ for \f$ x \notin V_b \f$.
 * */
template<typename scalar_type> scalar_type scalar_product(const vdom_vector<
		scalar_type>& v1, const vdom_vector<scalar_type>& v2) {
	if (v1.get_index_to_variable_id_map() == v2.get_index_to_variable_id_map()) {
		return scalar_product(v1.get_vector(),v2.get_vector());
	} else {
		index_to_variable_id_map_ptr new_iimap;
		index_type new_dim;
		index_to_index_bimap remap2;
		get_common_map(v1.get_index_to_variable_id_map(),
				v2.get_index_to_variable_id_map(), new_iimap, new_dim, remap2);
		/* only need to sum up variables that are both in l1 and l2. */
		scalar_type res = scalar_type(0);
		for (unsigned int i = 0; i < v2.size(); ++i) {
			unsigned int j = remap2.get_map(i);
			/* check whether it is in l1 and whether it is not zero */
			if (j < v1.size() && v1[j] != scalar_type(0)) {
				res += v1[j] * v2[i];
			}
		}
		//std::cout << "scalar multiply: " << l1 << " * " << l2 << " = " << res << std::endl;	
		return res;
	}
}
;

/** Product of vector elements.
 *
 * For two vectors \f$ a,b \f$ returns the vector \f$ c \f$  defined by \f$ c_x = a_x b_x \f$, if
 * \f$ x \in V_a \cap V_b \f$ and \f$ c_x=0 \f$ otherwise.
 * */
template<typename scalar_type> vdom_vector<scalar_type> element_product(
		const vdom_vector<scalar_type>& v1, const vdom_vector<scalar_type>& v2) {
	if (v1.get_index_to_variable_id_map() == v2.get_index_to_variable_id_map()) {
		return vdom_vector<scalar_type> (v1.domain(), element_product(
				v1.get_vector(), v2.get_vector()));
	} else {
		index_to_variable_id_map_ptr new_iimap;
		index_type new_dim;
		index_to_index_bimap remap2;
		get_common_map(v1.get_index_to_variable_id_map(),
				v2.get_index_to_variable_id_map(), new_iimap, new_dim, remap2);
		vdom_vector<scalar_type> res(new_iimap, scalar_type(0));
		/* only need to sum up variables that are both in l1 and l2. */
		for (unsigned int i = 0; i < v2.size(); ++i) {
			unsigned int j = remap2.get_map(i);
			/* check whether it is in l1 */
			if (j < v1.size()) {
				res[j] = v1[j] * v2[i];
			}
		}
		//std::cout << "element product: " << v1 << " * " << v2 << " = " << res << std::endl;
		return res;
	}
}
;

template<typename scalar_type> vdom_vector<scalar_type> operator*(
		const scalar_type& c, const vdom_vector<scalar_type>& v) {
	return vdom_vector<scalar_type> (c * v.get_vector(),
			v.get_index_to_variable_id_map());
}
;

template<typename scalar_type> vdom_vector<scalar_type> operator*(
		const vdom_vector<scalar_type>& v, const scalar_type& c) {
	return c * v;
}
;

template<typename scalar_type> vdom_vector<scalar_type> operator/(
		const vdom_vector<scalar_type>& v, const scalar_type& c) {
	return vdom_vector<scalar_type> (v.get_vector() / c,
			v.get_index_to_variable_id_map());
}
;

template<typename scalar_type> vdom_vector<scalar_type> operator-(
		const vdom_vector<scalar_type>& v) {
	return vdom_vector<scalar_type> (-v.get_vector(),
			v.get_index_to_variable_id_map());
}
;

template<typename scalar_type> vdom_vector<scalar_type>& operator+=(
		vdom_vector<scalar_type>& v1, const vdom_vector<scalar_type>& v2) {
	if (v1.get_index_to_variable_id_map() == v2.get_index_to_variable_id_map()) {
		v1.my_vector += v2.get_vector();
	} else {
		v1 = v1 + v2;
	}
	return v1;
}
;

template<typename scalar_type> vdom_vector<scalar_type>& operator*=(
		vdom_vector<scalar_type>& v, const scalar_type& c) {
	v.my_vector *= c;
	return v;
}
;

template<typename scalar_type> vdom_vector<scalar_type>& operator/=(
		vdom_vector<scalar_type>& v, const scalar_type& c) {
	v.my_vector /= c;
	return v;
}
;

template<typename scalar_type> vdom_vector<scalar_type>& operator/=(
		vdom_vector<scalar_type>& v, int c) {
	v /= scalar_type(c);
	return v;
}
;

//template<typename scalar_type>
inline vdom_vector<int>& operator/=(vdom_vector<int>& v, int c) {
	v /= c;
	return v;
}
;

/** Returns the elementwise absolute value of a vector. */
template<typename scalar_type> vdom_vector<scalar_type> vec_abs(const vdom_vector<
		scalar_type>& v) {
	return vdom_vector<scalar_type> (vec_abs(v.get_vector()),
			v.get_index_to_variable_id_map());
}
;

}

#endif /*vdom_VECTOR_OPERATORS_H_*/
