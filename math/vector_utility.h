#ifndef VECTOR_UTILITY_H_
#define VECTOR_UTILITY_H_

#include <iostream>
#include "math/vector.h"
#include "math/numeric/comp.h"
#include "type_conversion.h"
//#include "matrix.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "math/ublas_utility/ublas_utility.h"
#include "math/math_exception.h"

namespace math {

/** Returns a mapped vector.
 *
 * Returns a vector w defined by v[m(i)] for 0<=i<size().
 * Depending on m, w can be smaller or larger than v.
 * Elements not defined by the above are initialized with scalar_type(). */
template<typename scalar_type>
vector<scalar_type> map(const vector<scalar_type>& v,
		const index_to_index_bimap& m) {
	if (m.has_empty_codomain()) {
		return vector<scalar_type> ();
	} else {
		typename vector<scalar_type>::size_type m_i;
		vector<scalar_type> w(m.max_in_codomain() + 1);
		for (typename vector<scalar_type>::size_type i = 0; i < v.size(); ++i) {
			if (m.in_domain(i, m_i)) {
				w[m_i] = v[i];
			}
		}
		return w;
	}
}

/** Test if the 2D-points x,y,z are in counter-clockwise order.
 * The test is
 * x1 y2 + y1 z2 + x2 z1 - y2 z1 - x1 z2 - x2 y1 < 0.
 */
template<typename scalar_type> bool is_counter_clockwise_2D(const vector<
		scalar_type>& x, const vector<scalar_type>& y,
		const vector<scalar_type>& z) {
	assert(x.size()==2);
	assert(y.size()==2);
	assert(z.size()==2);

	// scalar_type d=x[0]*y[1]+y[0]*z[1]+x[1]*z[0]-y[1]*z[0]-x[0]*z[1]-x[1]*y[0];

	bool result = ((x[0] - y[0]) * (z[1] - y[1]) - (x[1] - y[1])
			* (z[0] - y[0]) < scalar_type(0));

	//std::cout << "ccw: " << x << "," << y << "," << z << ":" << result << std::endl;
	return result;

	//return d > scalar_type(0);
}
;

/** Test if the 3D-points x,y,z are in counter-clockwise order as seen from
 * the point p. */
template<typename scalar_type> bool is_counter_clockwise_3D(const vector<
		scalar_type>& x, const vector<scalar_type>& y,
		const vector<scalar_type>& z, const vector<scalar_type>& p) {
	assert(x.size()==3);
	assert(y.size()==3);
	assert(z.size()==3);
	assert(p.size()==3);

	using namespace boost::numeric::ublas;

	/** Test if det([x 1;y 1; z 1; p 1])>=0 */
	boost::numeric::ublas::matrix<scalar_type> M(4, 4);
	/*
	 project(M, range(0,1), range(0,3)) = x.get_vector_impl();
	 project(M, range(1,1), range(0,3)) = y.get_vector_impl();
	 project(M, range(2,1), range(0,3)) = z.get_vector_impl();
	 project(M, range(3,1), range(0,3)) = p.get_vector_impl();
	 */
	for (unsigned int i = 0; i < 3; ++i) {
		M(0, i) = x[i];
		M(1, i) = y[i];
		M(2, i) = z[i];
		M(3, i) = p[i];
	}
	for (unsigned int i = 0; i < 4; ++i) {
		M(i, 3) = scalar_type(1);
	}

	scalar_type d = boost::numeric::ublas::matrix_determinant(M);

	return d < scalar_type(0);
}
;

/** Returns the cross product of x and y (only defined in 3D). */
template<typename scalar_type> vector<scalar_type> cross_product_3D(
		const vector<scalar_type>& x, const vector<scalar_type>& y) {
	assert(x.size()==3);
	assert(y.size()==3);

	vector<scalar_type> z(3);
	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];
	return z;
}
;

/** Returns the angle between two vectors in 3D, with the vector ref defining
 *  the sign.
 *
 *  taken from http://www.gamedev.net/community/forums/topic.asp?topic_id=503639
 *  */
template<typename scalar_type> double angle_3D(const vector<scalar_type>& x,
		const vector<scalar_type>& y, const vector<scalar_type>& ref) {
	//				   vec3 c = cross(v1, v2);
	//				    float angle = std::atan2(length(c), dot(v1, v2));
	//				    return dot(c, reference) < 0.f ? -angle : angle;

	vector<scalar_type> c = cross_product_3D(x, y);
	double c_length = std::sqrt(convert_element<double> (scalar_product(c, c)));
	double dot = convert_element<double> (scalar_product(x, y));
	double s = convert_element<double> (scalar_product(c, ref));
	double angle = std::atan2(c_length, dot);
	if (s < 0.0)
		return -angle;
	else
		return angle;
}
;

/** A basic class for handling resolving exceptions. */
class resolving_exception: public math_exception {
public:
	resolving_exception(const std::string& msg) :
		math_exception(msg) {
	}
	;
};

/** Default conflict resolver for composing conflicting elements.
 *
 * Throws if the elements are not approximately equal, otherwise
 * keeps the first one. */
template<typename scalar_type> class throwing_resolver {
public:
	scalar_type resolve(const scalar_type& e1, const scalar_type& e2) const {
		if (!math::numeric::is_MEQ(e1, e2)) {
			std::string s = "Conflicting elements " + to_string(e1) + ","
					+ to_string(e2) + " found while composing.";
			throw resolving_exception(s);
			return scalar_type();
		}
		return e1;
	}
	;
};

/** Conflict resolver that memorizes whether a conflict has occurred.
 *
 * A conflict is when the elements are not approximately equal.
 * A conflicting element is set to scalar_type().
 *  */
template<typename scalar_type> class info_resolver {
public:
	info_resolver() : my_conflict(false) {};
	scalar_type resolve(const scalar_type& e1, const scalar_type& e2) const {
		if (!math::numeric::is_MEQ(e1, e2)) {
			info_resolver* nonconstthis = const_cast<info_resolver*>(this);
			nonconstthis->my_conflict = true;
			return scalar_type();
		}
		return e1;
	}
	;
	bool has_conflict() const {
		return my_conflict;
	}
private:
	bool my_conflict;
};

/** Map elements of vectors b1 and b2 onto a new vector using a conflict resolver
 * and return the resulting vector b.
 *
 * The conflict resolver determines what happens if the vectors have
 * an element that gets mapped to the same element in the result, and
 * these elements are different from each other.
 * The resolver is templated with scalar_type.
 * The default resolver throws if the values are not approximately
 * equal to each other.
 *
 * f1 and f2 map the elements of b1 and b2 to the resulting vector.
 *
 * If there is no conflict, then for all (i,j) in the index space of b1,
 *   b(f1(i))=b1(i),
 * and similarly for b2.
 */
template<typename scalar_type, template<typename > class resolver> vector<
		scalar_type> compose(const vector<scalar_type>& b1, const vector<
		scalar_type>& b2, const index_to_index_bimap& f1,
		const index_to_index_bimap& f2, const resolver<scalar_type>& resol = resolver<scalar_type>()) {
	typedef typename vector<scalar_type>::size_type size_type;
	size_type newdim=0;
	if(!f1.has_empty_codomain())
		newdim=std::max(newdim,f1.max_in_codomain() + 1);
	if(!f2.has_empty_codomain())
		newdim=std::max(newdim,f2.max_in_codomain() + 1);

	vector<scalar_type> b(newdim);

//	std::cout << "f1:" << f1 << ", f2:" << f2 << std::endl;

	/* Copy elements of b1 */
	for (size_type i = 0; i < b1.size(); ++i) {
		size_type f_i = f1.get_map(i);
		b[f_i] = b1[i];
	}

	/* Copy elements of b2 */
	for (size_type i = 0; i < b2.size(); ++i) {
		size_type f_i = f2.get_map(i);
		if (!f1.in_codomain(f_i))
			b[f_i] = b2[i];
		else {
			b[f_i] = resol.resolve(b[f_i], b2[i]);
		}
	}

	return b;
}
;

/** Map elements of vectors b1 and b2 onto a new vector b.
 *
 * Throws if the vectors have an element that gets mapped to the same
 * element in the result, and these elements are different from zero and from each other.
 * The resolver is templated with scalar_type.
 *
 * f1 and f2 map the elements of b1 and b2 to the resulting vector.
 *
 * If there is no conflict, then for all (i,j) in the index space of b1,
 *   b(f1(i))=b1(i),
 * and similarly for b2.
 */
template<typename scalar_type> vector<scalar_type> compose(const vector<
		scalar_type>& b1, const vector<scalar_type>& b2,
		const index_to_index_bimap& f1, const index_to_index_bimap& f2) {
	return compose<scalar_type, throwing_resolver> (b1, b2, f1, f2);
}
;

/** Returns the square of the 2-norm of a vector, which is
 * equal to the scalar product of the vector with itself.
 */
template<typename scalar_type> scalar_type norm2_square(const vector<
		scalar_type>& v) {
	scalar_type r(scalar_type(0));
	typedef typename vector<scalar_type>::size_type size_type;
	for (size_type i = 0; i < v.size(); ++i) {
		r += v[i] * v[i];
	}
	return r;
}
;

/** Returns the max of the elements in a vector.
 *
 * Throws if the vector is empty.
 */
template<typename scalar_type> scalar_type max(const vector<scalar_type>& v) {
	if (v.size() > 0) {
		scalar_type r = v[0];
		typedef typename vector<scalar_type>::size_type size_type;
		for (size_type i = 1; i < v.size(); ++i) {
			if (v[i] > r) {
				r = v[i];
			}
		}
		return r;
	} else
		return throw std::runtime_error("Cannot compute max of empty vector");
}
;

/** Returns the max of the elements in a vector, and assigns the corresponding index to k.
 *
 * Throws if the vector is empty.
 */
template<typename scalar_type> scalar_type max(const vector<scalar_type>& v,
		unsigned int& k) {
	if (v.size() > 0) {
		scalar_type r = v[0];
		k = 0;
		typedef typename vector<scalar_type>::size_type size_type;
		for (size_type i = 1; i < v.size(); ++i) {
			if (v[i] > r) {
				r = v[i];
				k = i;
			}
		}
		return r;
	} else {
		throw std::runtime_error("Cannot compute max of empty vector");
		return scalar_type(0);
	}
}
;

/** Returns the elementwise min of two vectors.
 *
 * Throws if the vectors are not of the same dimension.
 */
template<typename scalar_type> vector<scalar_type> min(const vector<scalar_type>& v1, const vector<scalar_type>& v2) {
	vector<scalar_type> res(v1);
	if (v1.size()==v2.size()) {
		typedef typename vector<scalar_type>::size_type size_type;
		for (size_type i = 0; i < v1.size(); ++i) {
			if (v2[i] < v1[i]) {
				res[i] = v2[i];
			}
		}
	} else
		throw std::runtime_error("Cannot compute min of vectors with different dimension");
	return res;
}
;

/** Returns the elementwise max of two vectors.
 *
 * Throws if the vectors are not of the same dimension.
 */
template<typename scalar_type> vector<scalar_type> max(const vector<scalar_type>& v1, const vector<scalar_type>& v2) {
	vector<scalar_type> res(v1);
	if (v1.size()==v2.size()) {
		typedef typename vector<scalar_type>::size_type size_type;
		for (size_type i = 0; i < v1.size(); ++i) {
			if (v2[i] > v1[i]) {
				res[i] = v2[i];
			}
		}
	} else
		throw std::runtime_error("Cannot compute max of vectors with different dimension");
	return res;
}
;

}

#endif /*VECTOR_UTILITY_H_*/
