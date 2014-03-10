#ifndef MATRIX_H_
#define MATRIX_H_

#include "global/global_types.h" // for __float128 abs

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include "math/basic_functions.h"
#include "math/vector.h"
#include "math/numeric/approx_comparator.h"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/storage.hpp>

namespace math {

using ::abs;

/** Matrix template class.
 * Based on class boost::numeric::ublas::matrix<scalar_type>
 * Requirements for scalar_type:
 * - needs to have a constructor that accepts scalar_type(0).
 * - needs to have a method abs().
 *  */

template<typename scalar_type> class matrix {
public:
	typedef boost::numeric::ublas::matrix<scalar_type> matrix_impl_type;
//	typedef boost::numeric::ublas::compressed_matrix<scalar_type> matrix_impl_type;
	typedef boost::numeric::ublas::matrix_row<matrix_impl_type> row_proxy_type;
	typedef boost::numeric::ublas::matrix_row<const matrix_impl_type> const_row_proxy_type;

	typedef scalar_type element_type;

	typedef typename matrix_impl_type::size_type size_type;
	typedef typename matrix_impl_type::difference_type difference_type;
	typedef typename matrix_impl_type::value_type value_type;
	typedef typename matrix_impl_type::const_reference const_reference;
	typedef typename matrix_impl_type::reference reference;
	//typedef typename matrix_impl_type::array_type array_type;
	typedef typename matrix_impl_type::const_closure_type const_closure_type;
	typedef typename matrix_impl_type::closure_type closure_type;
	typedef typename matrix_impl_type::vector_temporary_type vector_temporary_type;
	typedef typename matrix_impl_type::matrix_temporary_type matrix_temporary_type;
	typedef typename matrix_impl_type::storage_category storage_category;

	typedef typename matrix_impl_type::orientation_category orientation_category;

	/**
	 * Constructor list
	 *  */
	matrix() :
		my_matrix() {
	}

	matrix(size_type size1, size_type size2) :
		my_matrix(size1, size2, scalar_type(0)) {
	}

	/** Create a matrix with dimensions size1 by size2 and initialize all elements to default_value. */
	matrix(size_type size1, size_type size2, const scalar_type& default_value) :
		my_matrix(size1, size2) {

		for (unsigned int i = 0; i < my_matrix.size1(); i++) {
			for (unsigned int j = 0; j < my_matrix.size2(); j++)
				my_matrix(i, j) = default_value;
		}
	}

	explicit matrix(const matrix_impl_type &m) :
		my_matrix(m) {
	}

	matrix(const matrix &m) {
		my_matrix = m.get_matrix_impl();
	}

	/**
	 * Accessor
	 *  \return size1 or size2 (numbers of rows and columns of the matrix)
	 *  */
	size_type size1() const {
		return my_matrix.size1();
	}

	size_type size2() const {
		return my_matrix.size2();
	}

	/**
	 * Resizing
	 *  \param size1 or size2 (numbers of rows and columns of the matrix)
	 *  */
	void resize(size_type size1, size_type size2, bool preserve = true) {
		my_matrix.resize(size1, size2, preserve);
	}

	/**
	 * Element access
	 *  With operator()(size_type i, size_type j).
	 * \return reference or const_reference, depending on the context of the call
	 *  */
	const_reference operator()(size_type i, size_type j) const {
		return my_matrix(i, j);
	}

	reference operator()(size_type i, size_type j) {
		return my_matrix(i, j);
	}

	const_row_proxy_type row_proxy(size_type i) const {
		return const_row_proxy_type(my_matrix, i);
	}

	row_proxy_type row_proxy(size_type i) {
		return row_proxy_type(my_matrix, i);
	}

	vector<scalar_type> vector_from_row(size_type i) const {
		boost::numeric::ublas::matrix_row<const matrix_impl_type> mr(my_matrix,i);
		return vector<scalar_type>(mr);
	}

	vector<scalar_type> vector_from_column(size_type j) const {
		boost::numeric::ublas::matrix_column<const matrix_impl_type> mr(my_matrix,j);
//		vector<scalar_type> v(mr);
		vector<scalar_type> v(size1());
		for (size_t i=0;i<size1();++i) v[i]=my_matrix(i,j);
		assert(v.size()==size1());
		return v;
	}

	/**
	 *  Element assignment
	 * \param two size_type (information for inclusion) and a const_reference
	 * \return reference of element inserted
	 *  */
	reference insert_element(size_type i, size_type j, const_reference t) {
		return my_matrix.insert_element(i, j, t);
	}

	/**
	 *  Ecrase element
	 * \param two size_type (details of the element to ecrase)
	 *  */
	void erase_element(size_type i, size_type j) {
		my_matrix.erase_element(i, j);
	}

	/**
	 *  Zeroing matrix
	 *  */
	void clear() {
		my_matrix.clear();
	}

	/**
	 *  Matrix assignment
	 *  With operator= or function assign_temporary(matrix &m)
	 * \param Another matrix
	 *  */
	matrix& operator=(const matrix &m) {
		my_matrix = m.get_matrix_impl();
		return *this;
	}

	matrix& assign(const matrix &m) {
		my_matrix = m.get_matrix_impl();
		return *this;
	}

	/** Assign v to row i of *this. */
	void assign_row(size_type i,const vector<scalar_type>& v) {
		assert(i<size1());
		assert(size2()==v.size());
		row_proxy_type rowp(my_matrix,i);
		rowp=v.get_vector_impl();
	}

	/** Assign v to row i of *this. */
	void assign_row(size_type i,const const_row_proxy_type& v) {
		assert(i<size1());
		row_proxy_type rowp(my_matrix,i);
		rowp=v;
	}

	/** Assign v to column j of *this. */
	void assign_column(size_type j,const vector<scalar_type>& v) {
		assert(j<size2());
		assert(size1()==v.size());
		for(unsigned int i=0;i<size1();i++) {
			my_matrix(i,j) = v[i];
		}
	}

	/** Add a row to the end of the matrix. */
	void append_row(const vector<scalar_type>& x) {
		assert(size1()==0 || x.size()==size2());

		unsigned int new_row = size1();
		if (size1()==0) {
			my_matrix = matrix_impl_type(1, x.size());
		} else {
			unsigned int new_size = size1() + 1;
			my_matrix.resize(new_size, size2());
		}
		assign_row(new_row, x);
	}

	/** Add a matrix at the end of the matrix. */
	void append(const matrix<scalar_type>& M) {
		assert(size1()==0 || M.size2()==size2());

		unsigned int new_row = size1();
		if (size1()==0) {
			my_matrix = M.my_matrix;
		} else {
			unsigned int old_size = size1();
			unsigned int new_size = size1() + M.size1();
			my_matrix.resize(new_size, size2());
			submatrix_assign(M, old_size, new_size, 0, size2());
		}
	}

	matrix& assign_temporary(matrix &m) {
		swap(m);
		return *this;
	}

	/**
	 *  Matrix Swapping
	 *  */
	void swap(matrix &dm) {
		my_matrix.swap(dm.my_matrix);
	}

	friend void swap(matrix &m1, matrix &m2) {
		m1.swap(m2);
	}

	// Iterator types
	typedef typename matrix_impl_type::iterator1 iterator1;
	typedef typename matrix_impl_type::iterator2 iterator2;
	typedef typename matrix_impl_type::const_iterator1 const_iterator1;
	typedef typename matrix_impl_type::const_iterator2 const_iterator2;
	typedef typename matrix_impl_type::const_reverse_iterator1 const_reverse_iterator1;
	typedef typename matrix_impl_type::reverse_iterator1 reverse_iterator1;
	typedef typename matrix_impl_type::const_reverse_iterator2 const_reverse_iterator2;
	typedef typename matrix_impl_type::reverse_iterator2 reverse_iterator2;

	/**
	 *  Element lookup
	 *  begin, end, etc, and all find(int rank, size_type i, size_type j)
	 * \return iterator or const_iterator
	 *  */
	const_iterator1 find1(int rank, size_type i, size_type j) const {
		return my_matrix.find1(rank, i, j);
	}

	iterator1 find1(int rank, size_type i, size_type j) {
		return my_matrix.find1(rank, i, j);
	}

	const_iterator2 find2(int rank, size_type i, size_type j) const {
		return my_matrix.find2(rank, i, j);
	}

	iterator2 find2(int rank, size_type i, size_type j) {
		return my_matrix.find2(rank, i, j);
	}

	const_iterator1 begin1() const {
		return find1(0, 0, 0);
	}

	const_iterator1 end1() const {
		return find1(0, my_matrix.size1(), 0);
	}

	iterator1 begin1() {
		return find1(0, 0, 0);
	}

	iterator1 end1() {
		return find1(0, my_matrix.size1(), 0);
	}

	const_iterator2 begin2() const {
		return find2(0, 0, 0);
	}

	const_iterator2 end2() const {
		return find2(0, 0, my_matrix.size2());
	}

	iterator2 begin2() {
		return find2(0, 0, 0);
	}

	iterator2 end2() {
		return find2(0, 0, my_matrix.size2());
	}

	/**
	 *  ditto with reverse_iterator
	 *  */
	const_reverse_iterator1 rbegin1() const {
		return const_reverse_iterator1(end1());
	}

	const_reverse_iterator1 rend1() const {
		return const_reverse_iterator1(begin1());
	}

	reverse_iterator1 rbegin1() {
		return reverse_iterator1(end1());
	}

	reverse_iterator1 rend1() {
		return reverse_iterator1(begin1());
	}

	const_reverse_iterator2 rbegin2() const {
		return const_reverse_iterator2(end2());
	}

	const_reverse_iterator2 rend2() const {
		return const_reverse_iterator2(begin2());
	}

	reverse_iterator2 rbegin2() {
		return reverse_iterator2(end2());
	}

	reverse_iterator2 rend2() {
		return reverse_iterator2(begin2());
	}

	/**
	 * Get boost::numeric::ublas::matrix<scalar_type> member
	 * \return my_matrix
	 */
	const matrix_impl_type& get_matrix_impl() const {
		return my_matrix;
	}
	/**
	 * Gets the submatrix of A specified by the row and col bounds
	 * @return
	 */
	matrix<scalar_type> project_submatrix(unsigned int row_s, unsigned int row_e, unsigned int col_s, unsigned int col_e) const {
		namespace b = boost::numeric::ublas;
		b::range r(row_s,row_e);
		b::range c(col_s,col_e);
		math::matrix<scalar_type> my_new_matrix(b::project(my_matrix,r,c));

		return my_new_matrix;
	}
	/**
	 * Assign the submatrix with a new matrix
	 * @return
	 */
	void submatrix_assign(matrix<scalar_type> A, unsigned int row_s,unsigned int row_e, unsigned int col_s, unsigned int col_e){
		namespace b = boost::numeric::ublas;
		b::range r(row_s,row_e);
		b::range c(col_s,col_e);
		b::project(my_matrix,r,c) = A.get_matrix_impl();

	}
	/**
	 * \brief Computes infinity norm of a matrix.
	 *
	 * The infinity norm of a matrix is defined as the max row sum, taking absolute value of the row elements
	 * when summing. Returns zero if the matrix is empty.
	 * Infinity norm of a matrix A is defined as :
	 * ||A|| = max(\sum_{i=1}^cols |a_{i,j}|) over each row
	 */
	element_type infinity_norm() const {
		/*
		element_type max(0);
		if (size1() > 0 && size2() > 0) {

			unsigned int nrows = size1();
			unsigned int ncols = size2();

			element_type sum(0);
			for (unsigned int i = 0; i < nrows; i++) {
				sum = element_type(0);
				for (unsigned int j = 0; j < ncols; j++) {
					sum += abs_element(my_matrix(i, j));
				}
				if (sum > max) {
					max = sum;
				}
			}
			//u=norm_inf(my_matrix);
		}
		return max;
		*/
		return norm_inf(my_matrix);
	}
	;
	/**
	 * Returns the transpose matrix of *this
	 */
	matrix<element_type> transpose() const {
		/*
		 matrix<element_type> m(size2(), size1());
		 for (unsigned int i=0; i<m.size1(); i++) {
		 for (unsigned int j=0; j<m.size2(); j++)
		 m(i, j) = my_matrix(j,i);
		 }
		 return m;
		 */
		return matrix<element_type> (trans(my_matrix));
	}
	;
	/**
	 * \brief Computes the result of multiplication of *this matrix with element_type vector v
	 * and returns the result as a element_type vector.
	 *
	 * \param v element_type vector
	 * \return Multiplication result as a element_type vector
	 */
	vector<element_type> multiply_vector(const vector<element_type>& v) const {
		// Computes A * v
		if (size2() != v.size()) {
			std::stringstream ss;
			ss
					<< "Matrix-vector multiplication with wrong dimensions: Trying to multiply a "
					<< size1() << " by " << size2() << " matrix with a "
					<< v.size() << " by 1 vector.";
			throw std::runtime_error(ss.str());
		}
		/*
		 vector<element_type> result(size1());
		 element_type sum(0);
		 for (unsigned int i=0; i<size1(); i++) {
		 for (unsigned int j=0; j<size2(); j++)
		 sum = sum + my_matrix(i, j) * v[j];
		 result[i] = sum;
		 sum = element_type(0);
		 }
		 return result;
		 */
		return vector<element_type> (prod(my_matrix, v.get_vector_impl()));
	}
	;

	/**
	 * \brief Matrix scalar multiplication.
	 *
	 * \param delta element_type scalar value to be multiplied with each and every element
	 * of the matrix.
	 */
	void multiply_scalar(const element_type delta) {
		my_matrix *= delta;
		/*
		 for (unsigned int i=0; i<size1(); i++) {
		 for (unsigned int j=0; j<size2(); j++) {
		 my_matrix(i, j) = delta * my_matrix(i,j);
		 }
		 }
		 */
	}
	;
	/**
	 * \brief Multiplies a row of the matrix with a scalar
	 *
	 * \param delta Scalar value to be multiplied with the row elements
	 * \param i Row to be multiplied with the scalar.
	 */
	void multiply_row_scalar(const element_type delta, unsigned int i){
		for(unsigned int j=0;j<size2();j++)
			my_matrix(i,j) = delta * my_matrix(i,j);
	}

	/**
	 * \brief Multiplies a column of the matrix with a scalar
	 *
	 * \param delta Scalar value to be multiplied with the row elements
	 * \param j Column to be multiplied with the scalar.
	 */
	void multiply_col_scalar(const element_type delta, unsigned int j){
		for(unsigned int i=0;i<size1();i++)
			my_matrix(i,j) = delta * my_matrix(i,j);
	}

	/** Swaps rows i and j
	 *
	 * @attention Slow and inefficient implementation.
	 */
	void row_swap(size_t i, size_t j) {
		vector<scalar_type> old_row_i = vector_from_row(i);
		row_proxy_type rowi(my_matrix,i);
		rowi=row_proxy(j);
		assign_row(j,old_row_i);
	}

	/**
	 * Converts the elements of the vector to the type result_type.
	 * Uses the template function convert_element, which can be specialized to
	 * implement specific conversions between types.
	 */
	template<typename result_type> matrix<result_type> convert_to() const;

	/**
	 * \brief Checks if the matrix is a diagonal matrix.
	 *
	 * \return False if and only if the matrix has nonzero elements outside of its diagonal.
	 */
	bool is_diagonal() const {
		for (unsigned int i = 0; i < size1(); i++) {
			for (unsigned int j = 0; j < size2(); j++) {
				if (i != j) {
					if (!(my_matrix(i, j) == element_type(0)))
						return false;
				}
			}
		}
		return true;
	}

	/**
	 * Returns true if *this is a zero matrix.
	 *
	 * Returns true if all elements in the matrix are zero.
	 * An empty matrix is considered zero.
	 */
	bool is_zero() const {
		for (unsigned int i = 0;i < size1(); i++){
			for(unsigned int j = 0; j < size2(); j++){
					if(!(my_matrix(i,j) == element_type(0)))
						return false;
			}
		}
		return true;
	}

	/**
	 * Computes the inverse of the matrix.
	 *
	 * \return The inverse matrix if it exists. If the matrix is singular, then the parameter singular is set to true.
	 * singular is set to false otherwise.
	 * @notice This function works doesn't work for Rational type
	 */
	matrix<element_type> inverse(bool& singular) const {
		if (my_matrix.size1() != my_matrix.size2())
			std::runtime_error(
					"math::matrix:inverse()-> Inverse of non-square matrices not defined");

		if (my_matrix.size1() == 3 && my_matrix.size2() == 3) {
			const matrix_impl_type& A = my_matrix;
			element_type det = A(0, 0)
					* (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) - A(0, 1) * (A(1,
					0) * A(2, 2) - A(1, 2) * A(2, 0)) + A(0, 2) * (A(1, 0) * A(
					2, 1) - A(1, 1) * A(2, 0));
			singular = numeric::approx_comparator<element_type>::is_maybe_zero(det);
			if (singular) {
				return matrix<element_type>();
			} else {
				matrix<element_type> M(3, 3);
				matrix_impl_type& a = M.my_matrix;
				element_type s = element_type(1) / det;
				a(0, 0) = (s) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1));
				a(1, 0) = (s) * (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2));
				a(2, 0) = (s) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
				a(0, 1) = (s) * (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2));
				a(1, 1) = (s) * (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0));
				a(2, 1) = (s) * (A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1));
				a(0, 2) = (s) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
				a(1, 2) = (s) * (A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2));
				a(2, 2) = (s) * (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0));
				return M;
			}
		} else {
			using namespace boost::numeric::ublas;
			typedef permutation_matrix<std::size_t> pmatrix;

			// create a working copy of the input
			matrix_impl_type A(my_matrix);
			// create a permutation matrix for the LU-factorization
			pmatrix pm(A.size1());
			// perform LU-factorization
			int ret = lu_factorize(A, pm);
			if (ret != 0) {
				singular = true;
				return *this;
			}

			// create identity matrix of "inverse"
			matrix_impl_type temp(identity_matrix<element_type> (A.size1()));
			matrix_impl_type inverse(temp);

			// back-substitute to get the inverse
			lu_substitute(A, pm, inverse);
			singular = false;
			matrix<element_type> res(inverse);
			return res;
		}
	};

	/** Element-wise absolute
	 */
	matrix<element_type> abs() const {
		// create a working copy of the input
		matrix_impl_type A(size1(), size2());
		for (unsigned int i = 0; i < size1(); ++i) {
			for (unsigned int j = 0; j < size2(); ++j) {
				if (my_matrix(i, j)<element_type(0))
					A(i, j) = -(my_matrix(i, j));
				else
					A(i, j) = (my_matrix(i, j));
			}
		}
		matrix<element_type> res(A);
		return res;
	}
	;

	class output_format {
	public:
		output_format() { // default values
			preamble = "[\n";
			epilogue = "\n]";
			column_separator = ",";
			row_separator = "\n";
		}
		;
		static output_format matlab() {
			output_format f;
			f.preamble = "[";
			f.epilogue = "]";
			f.column_separator = ",";
			f.row_separator = ";";
			return f;
		}
		;
		std::string preamble; // written before matrix
		std::string epilogue; // written after matrix
		std::string column_separator; // written between columns (not after last)
		std::string row_separator; // written between rows (not after last)
	};
	static const output_format& get_output_format() {
		return my_output_format;
	}
	;
	static void set_output_format(output_format new_format) {
		my_output_format = new_format;
	}
	;

protected:
	matrix_impl_type my_matrix;
	static output_format my_output_format;
};

template<typename s> typename matrix<s>::output_format matrix<s>::my_output_format;

template<typename scalar_type, typename result_type> class matrix_converter {
public:
	static matrix<result_type> convert(const matrix<scalar_type>& m) {
		matrix<result_type> res(m.size1(), m.size2());
		for (unsigned int i = 0; i < m.size1(); i++) {
			for (unsigned int j = 0; j < m.size2(); j++) {
				res(i, j) = convert_element<result_type, scalar_type> (m(i, j));
			}
		}
		return res;
	}
	;
};

/** Partial specialization to avoid unnecessary casting if result_type is the
 * same as scalar_type. */
template<typename scalar_type> class matrix_converter<scalar_type, scalar_type> {
public:
	static matrix<scalar_type> convert(const matrix<scalar_type>& m) {
		return m;
	}
	;
};

template<typename scalar_type>
template<typename result_type> matrix<result_type> matrix<scalar_type>::convert_to() const {
	return matrix_converter<scalar_type, result_type>::convert(*this);
}
;

/*
 * Inverse method for Rational type matrices using naive LU decomposition
 */
template<>
inline matrix<Rational> matrix<Rational>::inverse(bool& singular) const {
	//Make a working copy of my_matrix
	std::size_t n = my_matrix.size1();
	std::size_t m = my_matrix.size2();
	if(n != m)
		throw std::runtime_error("Inverse not defined for non-square matrices");

	matrix<Rational> a(my_matrix);
	matrix<Rational> inverse(n,m);

	unsigned int k_prime;
	boost::numeric::ublas::vector<unsigned int> pi(n); // The permutation matrix

	for(unsigned int i = 0;i<n;i++)
		pi[i] = i;

	for(unsigned int k=0;k<n;k++){
		Rational p(0);
		for(unsigned int i=k; i<n; i++){
			if( abs_element(a(i, k)) > p){
				p = abs_element(a(i, k));
				k_prime = i;
			}
		}
		if(p == Rational(0))
		{
			singular = true;
			return *this;
		}
		unsigned int temp;
		temp = pi[k];
		pi[k] = pi[k_prime];
		pi[k_prime] = temp;

		for(unsigned int i =0;i<n;i++)
			a(k,i).swap(a(k_prime,i));
		for(unsigned int i=k+1;i<n;i++){
			a(i,k) = a(i,k)/a(k,k);
			for(unsigned int j=k+1;j<n;j++){
				a(i,j) = a(i,j) - a(i,k)*a(k,j);
			}
		}
	}

	// a is now decomposed into l and u matrix in place.

	/* LUx = Pb -> Ly = Pb
	 * We implement the LUP_Solve technique here, i.e, forward-backward substitution.
	 * We repeat the below process n times to get the inverse matrix.
	 */

	Rational sum;
	Rational b;
	boost::numeric::ublas::vector<Rational> y(n),x(n);
	for(unsigned int col = 0;col < m; col++){
		sum = Rational(0);
		for(int i=0;i<(int)n;i++){
			if(pi[i] == col)
				b=Rational(1); //??? check again
			else
				b=Rational(0);

			for(int j=0;j<=i-1;j++){
				sum = sum + a(i,j)*y(j);
			}
			y(i) = b - sum;
			sum = Rational(0);
		}
		// backward subst
		sum = Rational(0);
		for(int i = n-1;i>=0;i--){
			int j;
			for(j=i+1;j<(int)n;j++){
				sum = sum + a(i,j)*x(j);
			}
			inverse(i,col) = (y(i) - sum)/a(i,i);
			x(i) = inverse(i,col);
			sum = Rational(0);
		}
	}
	singular = false;
	return inverse;
}
;
}


/**
 * \brief Overloaded << operator for matrix class.
 */
template<typename scalar_type> std::ostream& operator<<(std::ostream& os,
		const math::matrix<scalar_type>& m) {
	os << m.get_output_format().preamble;
	for (unsigned int i = 0; i < m.size1(); i++) {
		if (i > 0) {
			os << m.get_output_format().row_separator;
		}
		for (unsigned int j = 0; j < m.size2(); j++) {
			if (j > 0) {
				os << m.get_output_format().column_separator;
			}
			os << m(i, j);
		}
	}
	os << m.get_output_format().epilogue;
	return os;
}
;

#include "math/matrix_operators.h"
#include "math/matrix_utility.h"

#endif /*MATRIX_H_*/
