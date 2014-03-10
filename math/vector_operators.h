#ifndef VECTOR_OPERATORS_H_
#define VECTOR_OPERATORS_H_

#include "math/numeric/approx_comparator.h"
#include "math/vector.h"

namespace math {

/** This file defines the basic vector operations.
 */

/** Element-wise comparison of two vectors. Considers two elements equal if
 * they are close, as defined by 
 * math::numeric::approx_comparator<scalar_type>::is_maybe_equal.
 * */
template<typename scalar_type> bool operator==(const vector<scalar_type>& v1,
		const vector<scalar_type>& v2) {
	if (v1.size() != v2.size())
		return false;
	typename vector<scalar_type>::const_iterator it1 = v1.begin();
	typename vector<scalar_type>::const_iterator it2 = v2.begin();
	while (it1 != v1.end()) {
		if (!math::numeric::approx_comparator<scalar_type>::is_maybe_equal(*it1,
				*it2))
			return false;
		++it1;
		++it2;
	}
	return true;
}
;

template<typename scalar_type> bool operator<(const vector<scalar_type>& v1,
		const vector<scalar_type>& v2) {
	if (v1.size() != v2.size())
		throw std::runtime_error("Comparing incomparable vectors");
	for (unsigned int i = 0; i < v1.size(); ++i) {
		if (!(v1[i] < v2[i]))
			return false;
		return true;
	}
	return true;
}
;

template<typename scalar_type> bool operator!=(const vector<scalar_type>& v1,
		const vector<scalar_type>& v2) {
	return !(v1 == v2);
}
;

template<typename scalar_type> vector<scalar_type> operator+(const vector<
		scalar_type>& v1, const vector<scalar_type>& v2) {
	assert(v1.size()==v2.size());
	/*
	 vector<scalar_type> res(v1.size());
	 typename vector<scalar_type>::const_iterator it=res.begin();
	 typename vector<scalar_type>::const_iterator it1=v1.begin();
	 typename vector<scalar_type>::const_iterator it2=v2.begin();
	 while (it1!=v1.end()) {
	 *it=*it1+*it2;
	 ++it;
	 ++it1;
	 ++it2;
	 }
	 */
	return vector<scalar_type> (v1.get_vector_impl() + v2.get_vector_impl());
}
;

template<typename scalar_type> vector<scalar_type> operator-(const vector<
		scalar_type>& v1, const vector<scalar_type>& v2) {
	assert(v1.size()==v2.size());
	return vector<scalar_type> (v1.get_vector_impl() - v2.get_vector_impl());
}
;

//template<typename scalar_type> scalar_type operator*(
//		const vector<scalar_type>& v1, const vector<scalar_type>& v2) {
//	assert(v1.size()==v2.size());
//	return boost::numeric::ublas::inner_prod(v1.get_vector_impl(),
//			v2.get_vector_impl());
//}
//;

/** Scalar product. */
template<typename scalar_type> scalar_type scalar_product(
		const vector<scalar_type>& v1, const vector<scalar_type>& v2) {
	assert(v1.size()==v2.size());
	return boost::numeric::ublas::inner_prod(v1.get_vector_impl(),
			v2.get_vector_impl());
}
;

/** Product of vector elements.
 *
 * For two vectors \f$ a,b \f$ of the same size returns the
 * vector \f$ c \f$  defined by \f$ c_i = a_i b_i \f$.
 * */
template<typename scalar_type> vector<scalar_type> element_product(
		const vector<scalar_type>& v1, const vector<scalar_type>& v2) {
	assert(v1.size()==v2.size());
	vector<scalar_type> res(v1);
	for (unsigned int i = 0; i < v1.size(); ++i) {
		res[i] *= v2[i];
	}
	//std::cout << "element multiply: " << v1 << " * " << v2 << " = " << res << std::endl;
	return res;

}
;

template<typename scalar_type> vector<scalar_type> operator*(
		const scalar_type& c, const vector<scalar_type>& v) {
	return vector<scalar_type> (c * v.get_vector_impl());
}
;

template<typename scalar_type> vector<scalar_type> operator*(const vector<
		scalar_type>& v, const scalar_type& c) {
	return c * v;
}
;

template<typename scalar_type> vector<scalar_type> operator/(const vector<
		scalar_type>& v, const scalar_type& c) {
	return vector<scalar_type> (v.get_vector_impl() / c);
}
;

template<typename scalar_type> vector<scalar_type> operator-(const vector<
		scalar_type>& v) {
	return vector<scalar_type> (-v.get_vector_impl());
}
;

template<typename scalar_type> vector<scalar_type>& operator+=(vector<
		scalar_type>& v1, const vector<scalar_type>& v2) {
	v1.my_vector += v2.get_vector_impl();
	return v1;
}
;

template<typename scalar_type> vector<scalar_type>& operator-=(vector<
		scalar_type>& v1, const vector<scalar_type>& v2) {
	v1.my_vector -= v2.get_vector_impl();
	return v1;
}
;

template<typename scalar_type> vector<scalar_type>& operator/=(vector<
		scalar_type>& v, const scalar_type& c) {
	v.my_vector /= c;
	return v;
}
;

template<typename scalar_type> vector<scalar_type>& operator*=(vector<
		scalar_type>& v, const scalar_type& c) {
	v.my_vector *= c;
	return v;
}
;

/** Returns the elementwise absolute value of a vector. */
template<typename scalar_type> vector<scalar_type> vec_abs(const vector<
		scalar_type>& v) {
	vector<scalar_type> res(v.size());
	for (unsigned int i = 0; i < v.size(); ++i) {
		if (v[i]>=scalar_type(0))
			res[i]=v[i];
		else
			res[i]=-v[i];
	}
	return res;
}
;



}

#endif /*VECTOR_OPERATORS_H_*/
