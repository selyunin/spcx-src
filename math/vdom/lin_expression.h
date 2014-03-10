#ifndef LIN_EXPRESSION_H_
#define LIN_EXPRESSION_H_

#include <boost/shared_ptr.hpp>

#include "math/vdom/vdom_vector.h"
#include "math/vdom/lin_expression_visitor.h"

namespace math {

/** A class for representing and computing with linear expressions of the
 * form
 *  a1 x1 + a2 x2 + ... + an xn + b.
 *  The variables x1,...,xn are associated to variable ids via an
 *  index_to_variable_id_map.
 *  The type of the coefficients a1,...,an,b is given by scalar_type.
 *
 * @attention The scalar_type must support the constructors scalar_type(0)
 * and scalar_type(1).
 *  */
template<typename scalar_type> class lin_expression: public printable {
public:
	typedef typename vdom_vector<scalar_type>::vector_type vector_type;
	typedef typename vector_type::size_type size_type;

	typedef typename vector_type::const_reference const_reference;
	typedef typename vector_type::reference reference;

	typedef typename vector_type::iterator iterator;
	typedef typename vector_type::const_iterator const_iterator;

	typedef typename vector_type::const_reverse_iterator const_reverse_iterator;
	typedef typename vector_type::reverse_iterator reverse_iterator;

	// --------------------------------------------
	/** \name Constructors
	 *  \{ */
	// --------------------------------------------

	/**
	 * Create a linear expression with no variables, and scalar b (by default 0).
	 *  */
	explicit lin_expression(scalar_type b = scalar_type(0)) :
		my_vec(), my_b(b) {
	}
	;

	/**
	 * Create a linear expression for a given vector of coefficients a1,...,an.
	 * The inhomogeneous coefficient b is set too scalar_type(0).
	 *  */
	explicit lin_expression(const vdom_vector<scalar_type> &v) :
		my_vec(v), my_b(0) {
	}
	;

	/**
	 * Create a linear expression for a given vector of coefficients a1,...,an,
	 * and a given scalar b.
	 *  */
	lin_expression(const vdom_vector<scalar_type> &v, const scalar_type& b) :
		my_vec(v), my_b(b) {
	}
	;

	/**
	 * Create a linear expression for a given vector of coefficients a1,...,an,
	 * a given scalar b and a given domain.
	 *  */
	lin_expression(const positional_vdomain& D, const vector_type &v, const scalar_type& b) :
		my_vec(D, v), my_b(b) {
	}
	;

	/**
	 * Create a linear expression for a given vector of coefficients a1,...,an,
	 * a given scalar b and a given index_to_variable_id_map.
	 *  */
	lin_expression(const vector_type &v, const scalar_type& b,
			const index_to_variable_id_map_ptr& pnew_map) :
		my_vec(v, pnew_map), my_b(b) {
	}
	;

	/**
	 * Create a linear expression for a given index_to_variable_id_map.
	 * The vector is allocated to the corresponding number of elements,
	 * but not initialized. The scalar b = 0.
	 *  */
	explicit lin_expression(const index_to_variable_id_map_ptr& pnew_map) :
		my_vec(pnew_map), my_b(0) {
	}
	;

	/** Copy constructor.
	 lin_expression(const lin_expression &v) :
	 my_vector(v.my_vector), my_b(v.my_b), my_iimap(v.my_iimap) {
	 }
	 ; */

	/** destructor. */
	~lin_expression() {
	}
	;

	/* \} */
	// --------------------------------------------
	/** \name Variable and domain methods
	 *  \{ */
	// --------------------------------------------

	/** Returns the \p variable_id at index i. */
	variable_id get_id(index_type i) const {
		return my_vec.get_id(i);
	}

	/** Returns the index of \p variable_id id. */
	bool has_id(variable_id id) const {
		return my_vec.has_id(id);
	}

	/** Returns the index of \p variable_id id. */
	index_type get_index(variable_id id) const {
		return my_vec.get_index(id);
	}

	/** Returns the index belonging to the id \p id if there is one.
	 * has_id is true if the id was found, and false otherwise.
	 * If has_id is false, then the returned index_type defaults to 0.
	 */
	index_type check_for_index(const variable_id& id, bool& has_id) const {
		return my_vec.check_for_index(id,has_id);
	}

	/** Returns the ids of all variables in the index_to_variable_id_map. */
	const variable_id_set& get_variable_ids() const {
		return my_vec.get_variable_ids();
	}

	/** Returns the ids of all variables with nonzero coefficients. */
	variable_id_set get_used_variable_ids() const {
		using namespace math::numeric;
		variable_id_set vis;
		for (size_t i = 0; i<my_vec.size(); ++i) {
			if (!is_MEQ(my_vec[i],scalar_type(0))) {
				vis.insert(get_id(i));
			}
		}
		return vis;
	}

	/** Returns the ids of the variables that are primed to degree \p
	 * prime_count. */
	variable_id_set get_primed_variables(unsigned int prime_count) const {
		return my_vec.get_primed_variables(prime_count);
	}

	/** Get the domain. */
	positional_vdomain domain() const {
		return my_vec.domain();
	}

	/** Get the variables in the domain. */
	std::set<variable> get_variables() const {
		return my_vec.get_variables();
	}

	/** Returns the position of variable x. */
	index_type pos(const variable& x) const {
		return my_vec.pos(x);
	}

	/** Returns the variable at position i. */
	variable get_variable(const index_type& i) const {
		return my_vec.get_variable(i);
	}

	/** Returns a pointer to the index_to_variable_id_map of \p *this. */
	const index_to_variable_id_map_ptr& get_index_to_variable_id_map() const {
		return my_vec.get_index_to_variable_id_map();
	}

	/* \} */
	// --------------------------------------------
	/** \name Domain modification
	 *  \{ */
	// --------------------------------------------

	/** Reorder *this according to the index_to_variable_id_map iimap.
	 * Throws if there is a nonzero coeffient for a variable
	 * that is not in iimap.
	 * */
	void reorder(const index_to_variable_id_map_ptr& iimap) {
		my_vec.reorder(iimap);
	}
	;

	/** Reorder *this according to the domain dom.
	 *
	 * Throws if there is a nonzero coeffients for a variable
	 * that is not in the domain. If dom has new variables, the
	 * corresponding coefficients are set to zero.
	 * */
	void reorder(const positional_vdomain& dom) {
		my_vec.reorder(dom);
	}
	;

	/** Remap *this according to the index_to_variable_id_map iimap.
	 * Any variable not in iimap is disregarded.
	 *
	 * Differs from reorder in that it doesn't throw. */
	void remap(const index_to_variable_id_map_ptr& iimap) {
		my_vec.remap(iimap);
	}
	;

	/** Remap *this according to the domain dom.
	 * Any variable not in dom is disregarded.
	 *
	 * Differs from reorder in that it doesn't throw. */
	void remap(const positional_vdomain& dom) {
		my_vec.remap(dom);
	}
	;

	/** Remove the variables in vis, i.e., remove them from the domain. */
	void remove_variables(const variable_id_set& vis) {
		my_vec.remove_variables(vis);
	}
	;

	/** Set the primedness of the variables with primedness of degree \p d to
	 * degree \p p. */
	void reassign_primedness(unsigned int d, unsigned int p = 0) {
		my_vec.reassign_primedness(d,p);
	}
	;

	/** Increase the primedness of the variables with primedness of
	 * degree \p d by 1.
	 * If d is 0, increase all. */
	void increase_primedness(unsigned int d = 0) {
		my_vec.reassign_primedness(d);
	}
	;

	/** Decrease the primedness of the variables with primedness of
	 * degree \p d by 1.
	 * If d is 0, decrease all. */
	void decrease_primedness(unsigned int d = 0) {
		my_vec.reassign_primedness(d);
	}
	;

	/* \} */
	// --------------------------------------------
	/** \name Element access
	 *  \{ */
	// --------------------------------------------

	/**
	 * Element access, e.g., v[i].
	 * \return reference or const_reference, depending on the context of the call
	 *  */
	const_reference operator [](size_type i) const {
		return my_vec[i];
	}
	;

	reference operator [](size_type i) {
		return my_vec[i];
	}
	;

	/**
	 * Element access, e.g., v(x), where x is a variable.
	 * \return reference or const_reference, depending on the context of the call
	 *  */
	const_reference operator()(const variable& x) const {
		return my_vec(x);
	}
	;

	reference operator()(const variable& x) {
		return my_vec(x);
	}
	;

	/** Returns the vdom_vector storing the homogenuous coefficients.
	 */
	const vdom_vector<scalar_type>& get_vdom_vec() const {
		return my_vec;
	}

	/** Returns the vdom_vector storing the homogenuous coefficients.
	 */
	vdom_vector<scalar_type>& get_vdom_vec() {
		return my_vec;
	}

	/**
	 * Get coefficients a1,...,an as vector<scalar_type>.
	 */
	const vector_type& get_vector() const {
		return my_vec.get_vector();
	}
	;

	/**
	 * Get the inhomogeneous coefficient b.
	 */
	const scalar_type& get_inh_coeff() const {
		return my_b;
	}
	;

	/**
	 * Set b.
	 */
	void set_inh_coeff(const scalar_type& b) {
		my_b = b;
	}
	;

	/** Get the coefficient of variable id to the value x.
	 * If id is not in the domain of *this, return zero.
	 */
	const scalar_type& get_coeff_with_id(const variable_id& id) const {
		return my_vec.get_coeff_with_id(id);
	}
	;

	/** Set the coefficient of variable x to the value s.
	 * Note that this might be slow due to internal reassignment of the iimap.
	 */
	void set_coeff(const variable& x, const scalar_type& s) {
		my_vec.set_coeff_with_id(x.get_id(),s);
	}
	;

	/** Set the coefficient of variable id to the value s.
	 * Note that this might be slow due to internal reassignment of the iimap.
	 */
	void set_coeff_with_id(const variable_id& id, const scalar_type& s) {
		my_vec.set_coeff_with_id(id,s);
	}
	;

	/** Set the coefficient of variable id to the value x.
	 * Throws if the id has no previous value (is not yet in the iimap). */
	void set_existing_coeff(const variable& x, const scalar_type& s) {
		my_vec.set_existing_coeff_with_id(x.get_id(),s);
	}
	;

	/** Set the coefficient of variable id to the value x.
	 * Throws if the id has no previous value (is not yet in the iimap). */
	void set_existing_coeff_with_id(const variable_id& id, const scalar_type& s) {
		my_vec.set_existing_coeff_with_id(id,s);
	}
	;

	/* \} */
	// --------------------------------------------
	/** \name Element information
	 *  \{ */
	// --------------------------------------------

	/** Returns the number of elements in the vector. */
	size_type size() const {
		return my_vec.size();
	}
	;

	/** Returns true iff all homogeneous coefficents a1,...,an are zero. */
	bool is_homogeneous_coeffs_zero() const {
		return my_vec.is_zero();
	}
	;

	/** Returns true iff all elements of *this are zero. */
	bool is_zero() const {
		return is_homogeneous_coeffs_zero() && my_b == scalar_type(0);
	}
	;

	/* \} */
	// --------------------------------------------
	/** \name Assignment
	 *  \{ */
	// --------------------------------------------

	/**
	 *  Zeroing all coefficients and b.
	 *  */
	void clear() {
		my_vec.clear();
		my_b = scalar_type(0);
	}
	;

	/**
	 *  Vector assignment
	 *  With operator= or function assign_temporary(vector &m)
	 * \param v vector
	 *  */
	lin_expression& operator =(const lin_expression& v) {
		my_vec = v.my_vec;
		my_b = v.my_b;
		return *this;
	}
	;

	lin_expression& assign_temporary(lin_expression& v) {
		this->swap(v);
		return *this;
	}
	;

	/**
	 *  Vector Swapping
	 *  */
	void swap(lin_expression& v) {
		if (this != &v) {
			my_vec.swap(v.my_vec);
			std::swap(my_b, v.my_b);
		}
	}
	;

	friend void swap(lin_expression &v1, lin_expression &v2) {
		v1.swap(v2);
	}
	;

	/** Process the visitor v. */
	void accept(lin_expression_visitor<scalar_type>& v) const {
		v.visit(*this);
	}
	;

	/**
	 * Converts the elements of the vector to the type result_type.
	 * Uses the template function convert_element, which can be specialized to
	 * implement specific conversions between types.
	 */
	template<typename result_type> lin_expression<result_type> convert_to() const {
		result_type b = convert_element<result_type, scalar_type> (my_b);
		lin_expression<result_type> res(
				my_vec.template convert_to<result_type>(),
				b);
		return res;
	}
	;

	/** Prints the linear expression to os. The default format is
	 * a1 x1 + a2 x2 + ... + an xn + b
	 * Nonzero elements are omitted.
	 */
	void print(std::ostream& os) const;

	/** Prints the linear expression to os, without the inhomogeneous coefficient.
	 */
	void print_hom(std::ostream& os) const;

	class output_format;
	static const output_format& get_output_format() {
		return my_output_format;
	}
	;
	static void set_output_format(output_format new_format) {
		my_output_format = new_format;
	}
	;

	/* \} */

protected:
	vdom_vector<scalar_type> my_vec;
	scalar_type my_b;
	static output_format my_output_format;
};

/* need to declare static members */
template<typename scalar_type> typename lin_expression<scalar_type>::output_format
		lin_expression<scalar_type>::my_output_format;

}

template<typename scalar_type>
class converter<math::lin_expression<scalar_type>, scalar_type> {
public:
	static math::lin_expression<scalar_type> convert(const scalar_type& x) {
		return math::lin_expression<scalar_type>(x);
	}
	;
};

template<>
class converter<math::lin_expression<__float128>, __float128> {
public:
	static math::lin_expression<__float128> convert(const __float128& x) {
		return math::lin_expression<__float128>(x);
	}
	;
};

#include "math/vdom/lin_expression_output.h"
#include "math/vdom/lin_expression_operators.h"

#endif /*LIN_EXPRESSION_H_*/
