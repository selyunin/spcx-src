/*
 * container_comp.h
 *
 *  Created on: Mar 31, 2010
 *      Author: frehse
 */

#ifndef CONTAINER_COMP_H_
#define CONTAINER_COMP_H_

#include "comp.h"
#include "math/vector.h"
#include "math/matrix.h"
#include "math/vdom/vdom_vector.h"

namespace math {
namespace numeric {

/** Extensions of comparisons to container types. */

/** Comparison for iterators over scalar_type.
 *  Corresponds to x1==x2 in the exact case. */
template<typename iterator_type1, typename iterator_type2>
bool is_MEQ(const iterator_type1& ib, const iterator_type1& ie,
		const iterator_type2& jb, const iterator_type2& je) {
	iterator_type1 i = ib;
	iterator_type2 j = jb;
	typedef typename iterator_type1::value_type scalar_type;
	while (i != ie && j != je) {
		if (!is_MEQ(*i, *j))
			return false;
		++i;
		++j;
	}
	// return true if both are at the end, and false otherwise
	return (i == ie && j == je);
}
;

/** Comparison for STL-like containers over scalar_type.
 *  Corresponds to x1==x2 in the exact case. */
template<typename scalar_type, template<typename> class container_type >
bool is_MEQ(const container_type<scalar_type>& v, const container_type<scalar_type>& w) {
	return is_MEQ(v.begin(),v.end(),w.begin(),w.end());
}
;

/** Check if two vectors are equal.
 *
 * @author Matthias Althoff, 10 May 2010*/
template<typename scalar_type> bool is_MEQ(const vector<scalar_type>& vector1,
		const vector<scalar_type>& vector2) {

	//check if both vectors have equal length
	if (vector1.size() == vector2.size()) {
		//init counter i
		unsigned int i = 0;

		//it is important to check first if i<vector1.size()
		while (i < vector1.size() && is_MEQ(vector1[i], vector2[i])) {
			i++;
		}
		if (i == vector1.size()) {
			return true;
		}
	}
	return false;
}
;

/** Check if all vector elements are equal to a scalar.
 **/
template<typename scalar_type> bool is_MEQ(const vector<scalar_type>& vector1,
		const scalar_type& a) {

	//init counter i
	unsigned int i = 0;

	//it is important to check first if i<vector1.size()
	while (i < vector1.size() && is_MEQ(vector1[i], a)) {
		i++;
	}
	if (i == vector1.size()) {
		return true;
	}
	return false;
}
;

/** Check if all matrix elements are equal to a scalar.
 **/
template<typename scalar_type> bool is_MEQ(const matrix<scalar_type>& A,
		const scalar_type& a) {
	for (unsigned int i=0;i<A.size1();++i) {
		for (unsigned int j=0;j<A.size2();++j) {
			if (!is_MEQ(A(i,j),a))
				return false;
		}
	}
	return true;
}
;

template <typename scalar_type>
bool is_definitely_LT(const scalar_type& x, const scalar_type& y) {
	return definitely(is_LT(x,y));
}

/** Lexicographical vector comparison using numeric comparisons. */
template <typename scalar_type, template<typename> class container_type >
bool is_lex_LT (const container_type<scalar_type>& v1, const container_type<scalar_type>& v2) {
	/*
	const_iterator it1=v1.begin();
	const_iterator it2=v2.begin();
	while (it1 != v1.end() && it2 != v2.end()) {
//std::cout << *it1 << " vs " << *it2 << std::endl;
		if (bool(is_LT(*it1,*it2))) return true;
		else if (bool(is_LT(*it2,*it1))) return false;
		else {
			++it1;
			++it2;
		}
	}
	if (it1 == v1.end() && it2 != v2.end())
		return true;
	else
		return false;
	*/
	return std::lexicographical_compare(v1.begin(),v1.end(),v2.begin(),v2.end(),is_definitely_LT<scalar_type>);
}

/** Lexicographical vector comparison using numeric comparisons.
 */
template <typename scalar_type, template<typename> class container_type >
class lex_comp_less {
public:
	typedef typename container_type<scalar_type>::const_iterator const_iterator;

	bool operator() (const container_type<scalar_type>& v1, const container_type<scalar_type>& v2) const {
		return is_lex_LT(v1,v2);
	}
};

/** Lexicographical vdom_vector comparison using numeric comparisons.
 *
 * Specialization for vdom_vectors that includes checking for the domains.
 */
template <typename scalar_type>
class lex_comp_less<scalar_type,math::vdom_vector > {
public:
	typedef typename math::vdom_vector<scalar_type>::const_iterator const_iterator;

	bool operator() (const math::vdom_vector<scalar_type>& v1, const math::vdom_vector<scalar_type>& v2) const {
		if (v1.domain()==v2.domain())
			return is_lex_LT(v1,v2);
		else {
			math::vdom_vector<scalar_type> v1_dom2=v1;
			v1_dom2.remap(v2.domain());
			math::vdom_vector<scalar_type> v2_dom1=v2;
			v2_dom1.remap(v1.domain());

			return is_lex_LT(v1_dom2,v2) && is_lex_LT(v1,v2_dom1);
		}
	}
};

/** Specialization of lex_comp_less for matrices */
template <typename scalar_type>
class lex_comp_less<scalar_type,matrix> {
public:
	//typedef typename matrix<scalar_type>::const_iterator const_iterator;

	bool operator() (const matrix<scalar_type>& v1, const matrix<scalar_type>& v2) const {
		return std::lexicographical_compare(
				v1.get_matrix_impl().data().begin(),
				v1.get_matrix_impl().data().end(),
				v2.get_matrix_impl().data().begin(),
				v2.get_matrix_impl().data().end(),
				is_definitely_LT<scalar_type> );
	}
};

}
}

#endif /* CONTAINER_COMP_H_ */
