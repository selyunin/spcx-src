#ifndef UBLAS_UTILITY_H_
#define UBLAS_UTILITY_H_

#include <boost/numeric/ublas/lu.hpp> 
#include <boost/numeric/ublas/matrix_expression.hpp>

namespace boost {
namespace numeric {
namespace ublas {

/** Matrix determinant. Taken from http://www.nabble.com/How-to-compute-determinant--td19385759.html */
template<class matrix_T> typename matrix_T::value_type LU_determinant(
		ublas::matrix_expression<matrix_T> const& mat_r) {
	typename matrix_T::value_type det = typename matrix_T::value_type(1);

	matrix_T mLu(mat_r() );
	ublas::permutation_matrix<std::size_t> pivots(mat_r().size1() );

	int is_singular = lu_factorize(mLu, pivots);

	if (!is_singular) {
		for (std::size_t i=0; i < pivots.size(); ++i) {
			if (pivots(i) != i)
				det *= typename matrix_T::value_type(-1);

			det *= mLu(i, i);
		}
	} else
		det = typename matrix_T::value_type(0);

	return det;
}
;

template<class matrix_T> inline typename matrix_T::value_type matrix_determinant_2D(
		ublas::matrix_expression<matrix_T> const& mat_r) {
	assert(mat_r().size1()==mat_r().size2());
	assert(mat_r().size1()==2);
	return (mat_r()(0, 0)*mat_r()(1, 1))-(mat_r()(0, 1)*mat_r()(1, 0));
}
;

template<typename value_type> inline value_type matrix_determinant_3D(
		const value_type& a, const value_type& b, const value_type& c,
		const value_type& d, const value_type& e, const value_type& f,
		const value_type& g, const value_type& h, const value_type& i) {
	return (a*e*i)-(a*f*h)-(b*d*i)+(b*f*g)+(c*d*h)-(c*e*g);
}
;

template<class matrix_T> inline typename matrix_T::value_type matrix_determinant_3D(
		ublas::matrix_expression<matrix_T> const& m) {
	assert(m().size1()==m().size2());
	assert(m().size1()==3);
	return (m()(0, 0)*m()(1, 1)*m()(2, 2))-(m()(0, 0)*m()(1, 2)*m()(2, 1))-(m()(0, 1)*m()(1, 0)
			*m()(2, 2))+(m()(0, 1)*m()(1, 2)*m()(2, 0))+(m()(0, 2)*m()(1, 0)*m()(2, 1))-(m()(0,
			2)*m()(1, 1)*m()(2, 0));
}
;

template<class matrix_T> typename matrix_T::value_type matrix_determinant_4D(
		ublas::matrix_expression<matrix_T> const& m) {
	assert(m().size1()==m().size2());
	assert(m().size1()==4);
	/* a00 a01 a02 a03
	 * a10 a11 a12 a13
	 * a20 a21 a22 a23
	 * a30 a31 a32 a33
	 */

	return (m()(0, 0)*matrix_determinant_3D(m()(1, 1), m()(1, 2), m()(1, 3), m()(2, 1),
			m()(2, 2), m()(2, 3), m()(3, 1), m()(3, 2), m()(3, 3))) - (m()(0, 1)
			*matrix_determinant_3D(m()(1, 0), m()(1, 2), m()(1, 3), m()(2, 0), m()(2, 2),
					m()(2, 3), m()(3, 0), m()(3, 2), m()(3, 3))) + (m()(0, 2)
			*matrix_determinant_3D(m()(1, 0), m()(1, 1), m()(1, 3), m()(2, 0), m()(2, 1),
					m()(2, 3), m()(3, 0), m()(3, 1), m()(3, 3))) - (m()(0, 3)
			*matrix_determinant_3D(m()(1, 0), m()(1, 1), m()(1, 2), m()(2, 0), m()(2, 1),
					m()(2, 2), m()(3, 0), m()(3, 1), m()(3, 2)));
}
;

/** @todo Why is there no implementation for higher dimensions? */
template<class matrix_T> typename matrix_T::value_type matrix_determinant(
		ublas::matrix_expression<matrix_T> const& m) {
	assert(m().size1()==m().size2());
	if (m().size1()==1)
		return m()(0, 0);
	else if (m().size1()==2)
		return matrix_determinant_2D(m);
	else if (m().size1()==3)
		return matrix_determinant_3D(m);
	else if (m().size1()==4)
		return matrix_determinant_4D(m);
	else
		return LU_determinant(m);

	// implicit cast to double gives compiler error	return LU_determinant(mat_r); 
}
;

}
}
}

#endif /*UBLAS_UTILITY_H_*/
