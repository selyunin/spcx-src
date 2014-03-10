#ifndef VECTOR_H_
#define VECTOR_H_

#include "global/global_types.h" // for __float128 abs

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <limits>
#include "math/vdom/index_to_index_bimap.h"
#include "math/basic_functions.h"
#include "math/type_conversion.h"
#include "math/scalar_types/rational.h"

namespace math {

using ::abs;

/** Forward declaration of vector. */
template<typename scalar_type> class vector;

/** Forward declaration of friend template functions. */
template<typename scalar_type> vector<scalar_type>& operator+=(vector<
		scalar_type>& v1, const vector<scalar_type>& v2);
template<typename scalar_type> vector<scalar_type>& operator-=(vector<
		scalar_type>& v1, const vector<scalar_type>& v2);
template<typename scalar_type> vector<scalar_type>& operator*=(vector<
		scalar_type>& v, const scalar_type& c);
template<typename scalar_type> vector<scalar_type>& operator/=(vector<
		scalar_type>& v, const scalar_type& c);

/** Vector template class.
 * Based on class boost::numeric::ublas::vector<scalar_type>
 *  */
template<typename scalar_type> class vector {
public:
	typedef scalar_type value_type;
	typedef boost::numeric::ublas::vector<scalar_type> vector_impl_type;
	typedef typename vector_impl_type::size_type size_type;
	typedef typename vector_impl_type::pointer pointer;
	typedef typename vector_impl_type::const_pointer const_pointer;
	typedef typename vector_impl_type::const_reference const_reference;
	typedef typename vector_impl_type::reference reference;
	typedef typename vector_impl_type::array_type array_type;

	typedef typename vector_impl_type::iterator iterator;
	typedef typename vector_impl_type::const_iterator const_iterator;

	typedef typename vector_impl_type::const_reverse_iterator
			const_reverse_iterator;
	typedef typename vector_impl_type::reverse_iterator reverse_iterator;

	/**
	 * Constructor list
	 *  */
	vector() :
		my_vector() {
	}

	vector(const vector &v) :
		my_vector(v.get_vector_impl()) {
	}

	explicit vector(size_type size) :
		my_vector(size) {
	}

	/** Create a vector with dimension size and initialize all elements to default_value. */
	explicit vector(size_type size, const scalar_type& default_value) :
		my_vector(size) {
		for (unsigned int i = 0; i < my_vector.size(); i++) {
			my_vector(i) = default_value;
		}
	}

	explicit vector(const vector_impl_type &v) :
		my_vector(v) {
	}

	explicit vector(const std::vector<scalar_type>& v) :
		my_vector(v.size() ) {
		for (unsigned int i=0;i<v.size();++i) {
			my_vector[i]=v[i];
		}
	}

	/** Virtual destructor. */
	virtual ~vector() {
	}
	;

	/**
	 * size()
	 *  \return vector size
	 *  */
	size_type size() const {
		return my_vector.size();
	}

	/** Reallocates the vector so that it can hold n elements.
	 *
	 * Erases or appends elements in order to bring the vector to the prescribed
	 * size. Appended elements copies of value_type().
	 * When p == false then existing elements are not preserved and elements
	 * will not appended as normal. Instead the vector is in the same state as
	 * that after an equivalent sizing constructor. */
	void resize(size_type size, bool preserve = true) {
		my_vector.resize(size, preserve);
	}

	/**
	 * Element support
	 *  Find an element from its position (size_type i).
	 * \return pointer or const_pointer, depending on the context of the call
	 *  */
	pointer find_element(size_type i) {
		return my_vector.find_element(i);
	}

	const_pointer find_element(size_type i) const {
		return my_vector.find_element(i);
	}

	/**
	 * Element access
	 *  With operator[](size_type i).
	 * \return reference or const_reference, depending on the context of the call
	 *  */
	const_reference operator [](size_type i) const {
		return my_vector[i];
	}

	reference operator [](size_type i) {
		return my_vector[i];
	}

	/**
	 *  Element assignment
	 * \param a size_type (position for inclusion) and a const_reference
	 * \return reference of element inserted
	 *  */
	reference insert_element(size_type i, const_reference t) {
		return my_vector.insert_element(i, t);
	}

	/**
	 *  Ecrase element
	 * \param a size_type (position of the element to ecrase)
	 *  */
	void erase_element(size_type i) {
		my_vector.erase_element(i);
	}

	/**
	 *  Zeroing vector
	 *  */
	void clear() {
		my_vector.clear();
	}

	/**
	 *  Vector assignment
	 *  With operator= or function assign_temporary(vector &m)
	 * \param Another vector
	 *  */
	vector& operator =(const vector& v) {
		my_vector = v.get_vector_impl();
		return *this;
	}
	/**
	 * Scalar multiplication of a vector
	 * \param A scalar value.
	 */

	vector& operator*(scalar_type k) {
		/*
		 for (unsigned int i=0; i<size(); i++) {
		 my_vector(i) = k * my_vector(i);
		 }
		 */
		my_vector *= k;
		return *this;
	}

	vector& assign_temporary(vector& v) {
		swap(v);
		return *this;
	}

	/**
	 *  Vector Swapping
	 *  */
	void swap(vector& v) {
		if (this != &v) {
			my_vector.swap(v.my_vector);
		}
	}

	friend void swap(vector &v1, vector &v2) {
		v1.swap(v2);
	}

	/**
	 * Get boost::numeric::ublas::vector<scalar_type> member
	 * \return my_vector
	 */
	const vector_impl_type& get_vector_impl() const {
		return my_vector;
	}

	/**
	 *  Element lookup
	 *  begin, end, etc, and all find(size_type i)
	 * \return iterator or const_iterator
	 *  */
	const_iterator find(size_type i) const {
		return my_vector.find(i);
	}

	iterator find(size_type i) {
		return my_vector.find(i);
	}

	const_iterator begin() const {
		return find(0);
	}

	const_iterator end() const {
		return find(my_vector.size());
	}

	iterator begin() {
		return find(0);
	}

	iterator end() {
		return find(my_vector.size());
	}

	/**
	 *  ditto with reverse_iterator
	 *  */
	const_reverse_iterator rbegin() const {
		return const_reverse_iterator(end());
	}

	const_reverse_iterator rend() const {
		return const_reverse_iterator(begin());
	}

	reverse_iterator rbegin() {
		return reverse_iterator(end());
	}

	reverse_iterator rend() {
		return reverse_iterator(begin());
	}
	/**
	 * Self assign the elements of v from position s to e(excluding e).
	 *
	 * @param v The vector from which to copy
	 * @param s The start point
	 * @param e The end point
	 */
	void subvector_assign(const vector<scalar_type>& v, std::size_t s, std::size_t e)
	{
		namespace b = boost::numeric::ublas;
		b::range r(s,e);
		b::project(my_vector,r) = v.get_vector_impl();
	}
	/**
	 * Self project the elements from position s to e (excluding e).
	 *
	 * @param s The start point
	 * @param e The end point
	 */
	vector<scalar_type> subvector_project(std::size_t s, std::size_t e)
	{
		namespace b = boost::numeric::ublas;
		b::range r(s,e);
		math::vector<scalar_type> my_new_vector(b::project(my_vector,r));

		return my_new_vector;
	}
	/**
	 * Computes the infinity norm of the vector, i.e., the largest absolute value of any element in the vector.
	 * Returns zero if the vector is empty.
	 *
	 * \param u A scalar_type variable to which the norm is assigned.
	 */
	scalar_type infinity_norm() const {
		scalar_type u = scalar_type(0);
		if (size() > 0) {
			u = abs_element(my_vector(0));
			for (unsigned int i = 1; i < size(); i++) {
				if (abs_element(my_vector(i)) > u) {
					u = abs_element(my_vector(i));
				}
			}
			// for some reason boost doesn't see the definition of std::abs for rationals
			//scalar_type test=std::abs(scalar_type(0)); //even though this line compiles just fine
			//u=norm_inf(my_vector);
		}
		return u;
	}
	;

	/**
	 * Computes the 1-norm of the vector.
	 *
	 * Returns sum of the absolute values of the elements in the vector.
	 * Returns zero if the vector is empty.
	 *
	 * \param u A scalar_type variable to which the norm is assigned.
	 */
	scalar_type one_norm() const {
		scalar_type u = scalar_type(0);
		if (size() > 0) {
			for (unsigned int i = 0; i < size(); i++) {
				u += abs_element(my_vector(i));
			}
		}
		return u;
	}
	;

	/**
	 * Computes the squared 2-norm of the vector.
	 *
	 * Returns sum of the squares of the elements in the vector.
	 * Returns zero if the vector is empty.
	 *
	 * \param u A scalar_type variable to which the norm is assigned.
	 *
	 * @note
	 * The square root is not taken so that scalar_type can
	 * be a data type that doesn't support square root.
	 */
	scalar_type two_norm_squared() const {
		scalar_type u = scalar_type(0);
		if (size() > 0) {
			for (unsigned int i = 0; i < size(); i++) {
				u += my_vector(i)*my_vector(i);
			}
		}
		return u;
	}
	;

	/**
	 * Computes scalar product with vector v.
	 *
	 * @param v vector
	 * @return scalar_type value containing the scalar product of *this with v
	 */
	scalar_type scalar_product(vector<scalar_type> v) const {
		if (size() != v.size())
			throw std::runtime_error(
					"Cannot compute scalar product of 2 vectors of different dimensions");
		/*
		 scalar_type s(0);
		 for(unsigned int i=0;i<size();i++)
		 s = s + v[i]*my_vector(i);
		 return s;
		 */
		return inner_prod(my_vector, v.my_vector);
	}
	/**
	 * Converts the elements of the vector to the type result_type.
	 * Uses the template function convert_element, which can be specialized to
	 * implement specific conversions between types.
	 */
	template<typename result_type> vector<result_type> convert_to() const {
		vector<result_type> res(size());
		for (unsigned int i = 0; i < size(); i++) {
			res[i] = convert_element<result_type, scalar_type> (my_vector(i));
		}
		return res;
	}
	;

	/** Returns true iff all elements of *this are zero. */
	bool is_zero() const {
		for (typename vector_impl_type::const_iterator it = my_vector.begin(); it
				!= my_vector.end(); ++it) {
			if (!(*it == scalar_type(0)))
				return false;
		}
		return true;
	}
	;

	class output_format;
protected:
	friend vector<scalar_type>& operator+=<> (vector<scalar_type>& v1,
			const vector<scalar_type>& v2);
	friend vector<scalar_type>& operator-=<> (vector<scalar_type>& v1,
			const vector<scalar_type>& v2);
	friend vector<scalar_type>& operator*=<> (vector<scalar_type>& v,
			const scalar_type& c);
	friend vector<scalar_type>& operator/=<> (vector<scalar_type>& v,
			const scalar_type& c);

	vector_impl_type my_vector;

};

template<typename scalar_type> class vector<scalar_type>::output_format {
public:
	output_format() { // default values
		preamble = "[";
		epilogue = "]";
		element_separator = ",";
	}
	;
	static output_format matlab() {
		output_format f;
		f.preamble = "[";
		f.epilogue = "]";
		f.element_separator = ",";
		return f;
	}
	;
	static output_format space_separated() {
		output_format f;
		f.preamble = "";
		f.epilogue = "";
		f.element_separator = " ";
		return f;
	}
	;
	std::string preamble; // written before matrix
	std::string epilogue; // written after matrix
	std::string element_separator; // written between elements
};

typedef vector<Rational> rational_vector;
typedef vector<double> double_vector;

}

/** Prints the vector v in the format
 * [0:element0,1:element1,...,n-1:element(n-1)]
 */
template<typename scalar_type> std::ostream& operator<<(std::ostream& os,
		const math::vector<scalar_type>& v) {
	os << "[";
	for (unsigned int i = 0; i < v.size(); i++) {
		if (i > 0) {
			os << ",";
		}
		os << v[i];
	}
	os << "]";
	return os;
}

#include "math/vector_operators.h"
#include "math/vector_utility.h"

#endif /*VECTOR_H_*/
