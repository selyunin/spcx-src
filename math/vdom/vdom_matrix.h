#ifndef vdom_MATRIX_H_
#define vdom_MATRIX_H_

#include "math/matrix.h"
#include "vdom_vector.h"
#include "positional_vdomain.h"

namespace math {

/** Matrix template class.
 * Based on class boost::numeric::ublas::matrix<scalar_type>
 * Requirements for scalar_type:
 * - needs to have a constructor that accepts scalar_type(0).
 * - needs to have a method abs().
 *  */
template<typename scalar_type> class vdom_matrix {
public:
	typedef matrix<scalar_type> matrix_impl_type;
	typedef scalar_type element_type;

	typedef typename matrix_impl_type::size_type size_type;
	typedef typename matrix_impl_type::value_type value_type;
	typedef typename matrix_impl_type::const_reference const_reference;
	typedef typename matrix_impl_type::reference reference;
	/**
	 * Constructor list
	 *  */
	vdom_matrix() :
		my_matrix() {
	}

	/** Create a matrix with codomain codom and domain dom and
	 * initialize all elements to default_value (0 if omitted). */
	vdom_matrix(const positional_vdomain& codom,
			const positional_vdomain& dom,
			const scalar_type& default_value = scalar_type(0)) :
		my_matrix(codom.size(), dom.size(), default_value), my_dom(dom), my_codom(codom) {
	}

	/** Create a matrix with codomain codom and domain dom and matrix values M. */
	vdom_matrix(const positional_vdomain& codom,
			const positional_vdomain& dom, const matrix_impl_type& M) :
		my_matrix(M), my_dom(dom), my_codom(codom) {
		assert(M.size1()==codom.size());
		assert(M.size2()==dom.size());
	}

	/** Returns the size of the domain. */
	size_type size1() const {
		return my_matrix.size1();
	}

	/** Returns the size of the codomain. */
	size_type size2() const {
		return my_matrix.size2();
	}

	/** Get the domain. */
	const positional_vdomain& domain() const {
		return my_dom;
	}

	/** Get the codomain. */
	const positional_vdomain& codomain() const {
		return my_codom;
	}

	/** Returns the element (i,j). */
	const_reference operator()(size_type i, size_type j) const {
		return my_matrix(i, j);
	}

	reference operator()(size_type i, size_type j) {
		return my_matrix(i, j);
	}

	/** Returns the element associated with the variables x and y. */
	const_reference operator()(const variable& x, const variable& y) const {
		return my_matrix(my_codom.pos(x), my_dom.pos(y));
	}

	reference operator()(const variable& x, const variable& y) {
		return my_matrix(my_codom.pos(x), my_dom.pos(y));
	}

	/** Returns the row (i,...) as a vector. */
	vector<scalar_type> vector_from_row(size_type i) const {
		return my_matrix.vector_from_row(i);
	}

	/** Returns the row (i,...) as a vdom_vector. */
	vdom_vector<scalar_type> vdom_vector_from_row(size_type i) const {
		return vdom_vector<scalar_type>(my_dom,my_matrix.vector_from_row(i));
	}

	/** Swap */
	void swap(vdom_matrix &dm) {
		my_matrix.swap(dm.my_matrix);
		::swap(my_dom, dm.my_dom);
		::swap(my_codom, dm.my_codom);
	}

	friend void swap(vdom_matrix &m1, vdom_matrix &m2) {
		m1.swap(m2);
	}

	/** Assign v to the row that corresponds to variable x in the codomain of *this.
	 * Returns the position of the row. */
	size_type assign_row(const variable& x, const vdom_vector<scalar_type>& v) {
		positional_vdomain::size_type pos;
		if (my_codom.in_domain(x, pos)) {
			if (my_dom == v.domain()) {
				my_matrix.assign_row(my_codom.pos(x), v.get_vector());
			} else
				throw std::runtime_error(
						"vdom_matrix: cannot assign row with different domain.");
		} else {
			throw std::runtime_error(
					"vdom_matrix: assignment to row for variable "
							+ x.get_name() + ", which is not in the codomain.");
		}
		return pos;
	}

	/** Reorder the elements of the matrix to obtain a given codomain and domain.
	 *
	 *
	 * Throws if there is a nonzero coeffients for a variable
	 * that is not in the domain or codomain, 
	 * unless the coefficent concerns the same variable in both domains. 
	 * If dom or codom have new variables, the
	 * corresponding coefficients are set to zero.
	 * */
	void reorder(const positional_vdomain& codom,const positional_vdomain& dom) {
		if (domain() != dom || codomain() != codom) {
			matrix_impl_type M(codom.size(),dom.size());
			// f maps rows
			position_map f=compute_reordering(codomain(),codom,true);
			// g maps columns
			position_map g = compute_reordering(domain(), dom, true);
			for (positional_vdomain::size_type i = 0; i < codomain().size();
					i++) {
				positional_vdomain::size_type mapped_row = f.get_map(i);
				if (mapped_row != positional_vdomain::invalid_pos()) {
					for (positional_vdomain::size_type j = 0;
							j < domain().size(); j++) {
						positional_vdomain::size_type mapped_col = g.get_map(j);
						if (mapped_col != positional_vdomain::invalid_pos()) {
							M(mapped_row, mapped_col) = my_matrix(i, j);
						} else {
							// only throw if it's nonzero and doesn't affect the same variable
							if (my_matrix(i, j) != scalar_type(0) && codomain().get_variable(i)!=domain().get_variable(j)) {
								std::string varname = domain().get_variable(j).get_name();
								throw basic_exception(
										"Tried projecting away domain variable "+varname+" with nonzero coefficients.");
							}
						}
					}
				} else {
					// test if there are nonzero coefficients
					for (positional_vdomain::size_type j = 0;
							j < domain().size(); j++) {
						if (my_matrix(i, j) != scalar_type(0)) {
							std::string varname = domain().get_variable(j).get_name();
							throw basic_exception(
									"Tried projecting away codomain variable "+varname+" with nonzero coefficients.");
						}
					}
				}
			}
			vdom_matrix res(codom, dom, M);
			swap(res);
		}
	}
	;

	/**
	 * \brief Computes infinity norm of a matrix.
	 *
	 * The infinity norm of a matrix is defined as the max row sum, taking absolute value of the row elements
	 * when summing. Returns zero if the matrix is empty.
	 * Infinity norm of a matrix A is defined as :
	 * ||A|| = max(\sum_{i=1}^cols !a_{i,j}|) over each row
	 */
	element_type infinity_norm() const {
		return my_matrix.infinity_norm();
	}
	;
	/**
	 * Returns the transpose matrix of *this
	 */
	vdom_matrix<element_type> transpose() const {
		return vdom_matrix<element_type> (my_dom, my_codom,
				my_matrix.transpose());
	}
	;

	/* Element-wise absolute */
	vdom_matrix<element_type> abs() const {
		return vdom_matrix<element_type> (my_dom, my_codom,
				my_matrix.abs());
	}
	;
	/**
	 * \brief Computes the result of multiplication with vector v.
	 *
	 * Any coefficients for variables in the domain of v that are not
	 * in the domain of *this are ignored. If a variable from the domain
	 * of *this is not in the domain of v, v is considered to have a
	 * coefficient zero for that variable.
	 *
	 */
	vdom_vector<element_type> multiply_vector(
			const vdom_vector<element_type>& v) const {
		if (my_dom == v.domain()) {
			return vdom_vector<element_type> (my_codom,
					my_matrix.multiply_vector(v.get_vector()));
		} else {
			vdom_vector<element_type> result(my_codom, element_type(0));
			for (unsigned int j = 0; j < size2(); j++) {
				index_type j_mapped;
				if (v.domain().in_domain(my_dom.get_variable(j), j_mapped)) {
					for (unsigned int i = 0; i < size1(); i++) {
						result[i] += my_matrix(i, j) * v[j_mapped];
					}
				}
			}
			return result;
		}
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
	}
	;

	/**
	 * Converts the elements of the vector to the type result_type.
	 * Uses the template function convert_element, which can be specialized to
	 * implement specific conversions between types.
	 */
	template<typename result_type> vdom_matrix<result_type> convert_to() const;

	/**
	 * \brief Checks if the matrix is a diagonal matrix.
	 *
	 * \return False if and only if the matrix has nonzero elements outside of its diagonal.
	 */
	bool is_diagonal() const {
		return my_matrix.is_diagonal();
	}

	/**
	 * Returns true if *this is a zero matrix.
	 *
	 * Returns true if all elements in the matrix are zero.
	 * An empty matrix is considered zero.
	 */
	bool is_zero() const {
		return my_matrix.is_zero();
	}

	/*
	 * \brief Computes the inverse of the matrix.
	 *
	 * \return The inverse matrix if it exists. If the matrix is singular, then the parameter singular is set to true.
	 * singular is set to false otherwise.
	 */
	vdom_matrix<element_type> inverse(bool& singular) const {
		return vdom_matrix<element_type> (my_codom, my_dom, my_matrix.inverse(
				singular));
	}

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

	/**
	 * Get a normal matrix.
	 */
	const matrix_impl_type& get_matrix() const {
		return my_matrix;
	}
	;

	/**
	 * Get a normal matrix.
	 */
	matrix_impl_type& get_matrix() {
		return my_matrix;
	}
	;

protected:

	matrix_impl_type my_matrix;
	positional_vdomain my_dom;
	positional_vdomain my_codom;
	static output_format my_output_format;
};

template<typename s> typename vdom_matrix<s>::output_format
		vdom_matrix<s>::my_output_format;

template<typename scalar_type>
template<typename result_type> vdom_matrix<result_type> vdom_matrix<
		scalar_type>::convert_to() const {
	return vdom_matrix<result_type> (my_codom, my_dom, my_matrix.convert_to<
			result_type> ());
}
;

}

/**
 * \brief Overloaded << operator for matrix class.
 */
template<typename scalar_type> std::ostream& operator<<(std::ostream& os,
		const math::vdom_matrix<scalar_type>& m) {
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

#endif /*vdom_MATRIX_H_*/
