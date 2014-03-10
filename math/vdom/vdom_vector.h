#ifndef vdom_VECTOR_H_
#define vdom_VECTOR_H_

#include "math/vector.h"
#include <boost/shared_ptr.hpp>
#include "math/vdom/index_to_variable_id_map.h"
#include "math/vdom/variable.h"
#include "utility/printable.h"
#include "math/vdom/index_to_variable_id_map_provider.h"

namespace math {

using ::operator<<;

/** Forward declaration. */
template<typename scalar_type> class vdom_vector;

/** Forward declaration of friend template functions. */
template<typename scalar_type> vdom_vector<scalar_type>& operator+=(
		vdom_vector<scalar_type>& v1, const vdom_vector<scalar_type>& v2);
template<typename scalar_type> vdom_vector<scalar_type>& operator*=(
		vdom_vector<scalar_type>& v, const scalar_type& c);
template<typename scalar_type> vdom_vector<scalar_type>& operator/=(
		vdom_vector<scalar_type>& v, const scalar_type& c);

/** Visitor for processing vdom_vector objects */
template<typename scalar_type> class vdom_vector_visitor {
public:
	typedef unsigned int index_type;
	typedef scalar_type element_type;
	virtual ~vdom_vector_visitor() {
	}
	;
	/** visit is called by the vdom_vector's accept. */
	virtual void visit(const vdom_vector<scalar_type>& l) {
		vdom_vector_prologue(l);
		visit_coeffs(l);
		vdom_vector_epilogue(l);
	}
	;
	/** The prologue is called before processing the linear expression. */
	virtual void vdom_vector_prologue(const vdom_vector<scalar_type>& l) {
	}
	;
	/** If it returns true the visit_coeff function will be called for every coefficient. */
	virtual void visit_coeffs(const vdom_vector<scalar_type>& l) {
		for (unsigned int i = 0; i < l.size(); i++) {
			visit_coeff(l, l[i], i);
		}
	}
	;
	/** The visit_coeff is called for every coefficient with its index i and the iimap of the vdom_vector. */
	virtual void visit_coeff(const vdom_vector<scalar_type>& l,
			const scalar_type& coeff, const index_type& i) {
	}
	;
	/** The epilogue is called after processing the linear expression. */
	virtual void vdom_vector_epilogue(const vdom_vector<scalar_type>& l) {
	}
	;
};

/** A class for representing and computing with vectors whose values are
 *  associated with named variables.
 *  The variables x1,...,xn are associated to variable ids via an
 *  index_to_variable_id_map.
 *  The type of the coefficients (values) a1,...,an,b is given by scalar_type.
 *
 * @attention The scalar_type must support the constructors scalar_type(0)
 * and scalar_type(1).
 *  */
template<typename scalar_type> class vdom_vector: public index_to_variable_id_map_provider,
		public printable {
public:
	typedef scalar_type value_type;
	typedef math::vector<scalar_type> vector_type;
	typedef typename vector_type::size_type size_type;

	typedef typename vector_type::const_reference const_reference;
	typedef typename vector_type::reference reference;

	typedef typename vector_type::iterator iterator;
	typedef typename vector_type::const_iterator const_iterator;

	typedef typename vector_type::const_reverse_iterator const_reverse_iterator;
	typedef typename vector_type::reverse_iterator reverse_iterator;

	/**
	 * Create a linear expression with no variables, and scalar b=0.
	 *  */
	vdom_vector() :
		my_vector() {
	}
	;

	/**
	 * Create a zero vector for domain D. */
	vdom_vector(const positional_vdomain& D) : index_to_variable_id_map_provider(D.my_iimap),
		my_vector(D.size(),scalar_type(0)) {
		assert(D.my_iimap->dimensions()==D.size());
	}
	;

	/**
	 * Create a named vector for a given vector of coefficients a1,...,an,
	 * and a given domain. */
	vdom_vector(const positional_vdomain& D, const vector_type& v) : index_to_variable_id_map_provider(D.my_iimap),
		my_vector(v) {
		assert(D.my_iimap->dimensions()==v.size());
	}
	;

	/**
	 * Create a named vector for a given domain and initialize all elements
	 * to e. */
	vdom_vector(const positional_vdomain& D, const scalar_type& e) : index_to_variable_id_map_provider(D.my_iimap),
			my_vector(vector_type(D.my_iimap->dimensions(),e))	{
	}
	;

	/**
	 * Create a named vector for a given vector of coefficients a1,...,an,
	 * and a given index_to_variable_id_map. */
	vdom_vector(const vector_type &v,
			const index_to_variable_id_map_ptr& pnew_map) : index_to_variable_id_map_provider(pnew_map),
		my_vector(v) {
		assert(pnew_map->dimensions()==v.size());
	}
	;

	/**
	 * Create a named vector for a given index_to_variable_id_map.
	 * The vector is allocated to the corresponding number of elements,
	 * and is initialized to scalar_type(0).
	 *  */
	vdom_vector(const index_to_variable_id_map_ptr& pnew_map) : index_to_variable_id_map_provider(pnew_map),
		my_vector(pnew_map->dimensions(), scalar_type(0)) {
	}
	;

	/** Copy constructor.
	 *
	 * It's different from the default one, since it has to
	 * @remark call the base class copy constructor as well.
	 */
	vdom_vector(const vdom_vector& v) : index_to_variable_id_map_provider(v),
		my_vector(v.my_vector) {
		assert(v.domain().size()==v.my_vector.size());
	}
	;

	/** Virtual destructor. */
	virtual ~vdom_vector() {
	}
	;

	/** Returns the number of elements in the vector. */
	size_type size() const {
		return my_vector.size();
	}
	;

	/**
	 * Element access, e.g., v[i].
	 * \return reference or const_reference, depending on the context of the call
	 *  */
	const_reference operator [](size_type i) const {
		return my_vector[i];
	}
	;

	reference operator [](size_type i) {
		return my_vector[i];
	}
	;

	/**
	 * Element access, e.g., v(x), where x is a variable.
	 * \return reference or const_reference, depending on the context of the call
	 *  */
	const_reference operator()(const variable& x) const {
		return my_vector[domain().pos(x)];
	}
	;

	reference operator()(const variable& x) {
		return my_vector[domain().pos(x)];
	}
	;

	/** Returns true iff all elements of *this are zero. */
	virtual bool is_zero() const {
		return my_vector.is_zero();
	}
	;

	/**
	 *  Element assignment
	 * \param i size_type (position for inclusion)
	 * \param t const_reference
	 * \return reference of element inserted
	 *  */
	reference insert_element(size_type i, const_reference t) {
		return my_vector.insert_element(i, t);
	}
	;

	/**
	 *  Erase element
	 * \param i size_type (position of the element to ecrase)
	 *  */
	void erase_element(size_type i) {
		my_vector.erase_element(i);
	}
	;

	/**
	 *  Set all coefficients to zero, leaving variable placement intact.
	 *  */
	virtual void clear() {
		my_vector.clear();
	}
	;

	/**
	 *  Vector assignment
	 *  With operator= or function assign_temporary(vector &m)
	 * \param v vector
	 *  */
	vdom_vector& operator =(const vdom_vector& v) {
		index_to_variable_id_map_provider::operator=(v);
		my_vector = v.my_vector;
		return *this;
	}
	;

	vdom_vector& assign_temporary(vdom_vector& v) {
		swap(v);
		return *this;
	}
	;

	/**
	 *  Vector Swapping
	 *  */
	void swap(vdom_vector& v) {
		if (this != &v) {
			index_to_variable_id_map_provider::swap(v);
			my_vector.swap(v.my_vector);
		}
	}
	;

	friend void swap(vdom_vector &v1, vdom_vector &v2) {
		v1.swap(v2);
	}
	;

	/**
	 * Get coefficients a1,...,an as vector<scalar_type>.
	 */
	const vector_type& get_vector() const {
		return my_vector;
	}
	;

	/** Get the coefficient of variable id to the value x.
	 * If id is not in the domain of *this, return zero.
	 */
	const scalar_type& get_coeff_with_id(const variable_id& id) const {
		bool found;
		index_type i = get_index_to_variable_id_map()->check_for_index(id,
				found);
		if (!found) {
			return my_zero;
		} else {
			return my_vector[i];
		}
	}
	;

	/**
	 * Set coefficients a1,...,an as vector<scalar_type>.
	 */
	void set_vector(const vector_type& v) {
		my_vector = v;
	}
	;

	/** Set the coefficient of variable v to the value s.
	 * Note that this might be slow due to internal reassignment of positional_vdomain.
	 */
	void set_coeff_with_id(const variable& v, const scalar_type& s) {
		variable_id id = v.get_id();
		set_coeff_with_id(id,s);
	}
	;

	/** Set the coefficient of variable id to the value s.
	 * Note that this might be slow due to internal reassignment of the iimap.
	 */
	void set_coeff_with_id(const variable_id& id, const scalar_type& s) {
		bool found;
		index_type i = get_index_to_variable_id_map()->check_for_index(id,
				found);
		if (!found) {
			set_index_to_variable_id_map(
					get_index_to_variable_id_map()->get_map_with_id_added(id));
			i = get_index(id);
			my_vector.resize(i + 1);
		}
		my_vector[i] = s;
	}
	;

	/** Set the coefficient of variable id to the value x.
	 * Throws if the id has no previous value (is not yet in the iimap). */
	void set_existing_coeff_with_id(const variable_id& id, const scalar_type& s) {
		my_vector[get_index_to_variable_id_map()->get_index(id)] = s;
	}
	;

	/**
	 *  Element lookup
	 *  begin, end, etc, and all find(size_type i)
	 * \return iterator or const_iterator
	 *  */
	const_iterator find(size_type i) const {
		return my_vector.find(i);
	}
	;

	iterator find(size_type i) {
		return my_vector.find(i);
	}
	;

	const_iterator begin() const {
		return find(0);
	}
	;

	const_iterator end() const {
		return find(my_vector.size());
	}
	;

	iterator begin() {
		return find(0);
	}
	;

	iterator end() {
		return find(my_vector.size());

	}
	;

	/**
	 *  ditto with reverse_iterator
	 *  */
	const_reverse_iterator rbegin() const {
		return const_reverse_iterator(end());
	}
	;

	const_reverse_iterator rend() const {
		return const_reverse_iterator(begin());
	}
	;

	reverse_iterator rbegin() {
		return reverse_iterator(end());
	}
	;

	reverse_iterator rend() {
		return reverse_iterator(begin());
	}
	;

	/**
	 * Computes the infinity norm of the vector, i.e., the largest absolute value of any element in the vector.
	 * Returns zero if the vector is empty.
	 *
	 * \param u A scalar_type variable to which the norm is assigned.
	 */
	scalar_type infinity_norm() const {
		return my_vector.infinity_norm();
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
		return my_vector.one_norm();
	}
	;
	/**
	 * Computes the scalar product of the two vectors.
	 *
	 * @param v The vector to scalar_product
	 * @return The scalar product of this with v.
	 */
	scalar_type scalar_product(const vdom_vector& v) const {
		return my_vector.scalar_product(v.get_vector());
	}
	;
	/** Reorder *this according to the index_to_variable_id_map iimap.
	 * Throws if there is a nonzero coefficient for a variable
	 * that is not in iimap.
	 * @todo: optimize by swapping elements instead of creating a new vector
	 * @todo: ensure that all elements of *this are in the iimap
	 * */
	virtual void reorder(const index_to_variable_id_map_ptr& iimap) {
		if (iimap != get_index_to_variable_id_map()) {
			vdom_vector res(iimap);
			for (unsigned int i = 0; i < size(); i++) {
				if (!(my_vector[i] == scalar_type(0))) {
					res.set_existing_coeff_with_id(get_id(i), my_vector[i]);
				}
			}
			*this = res;
		}
	}
	;

	/** Reorder *this according to the domain dom.
	 *
	 * Throws if there is a coefficient different from default_value for a variable
	 * that is not in the domain. If dom has new variables, the
	 * corresponding coefficients are set to zero.
	 * */
	virtual void reorder(const positional_vdomain& dom, const scalar_type& default_value = scalar_type(0)) {
		const index_to_variable_id_map_ptr& iimap = dom.get_index_to_variable_id_map();
		if (iimap != get_index_to_variable_id_map()) {
			vdom_vector res(iimap, default_value);
			for (unsigned int i = 0; i < size(); i++) {
				if (!(my_vector[i] == default_value)) {
					res.set_existing_coeff_with_id(get_id(i), my_vector[i]);
				}
			}
			*this = res;
		}
	}
	;


	/** Remap *this according to the index_to_variable_id_map iimap.
	 * Any variable not in iimap is disregarded.
	 *
	 * Differs from reorder in that it doesn't throw. */
	virtual void remap(const index_to_variable_id_map_ptr& iimap) {
		if (iimap != get_index_to_variable_id_map()) {
			vdom_vector res(iimap);
			bool found;
			for (unsigned int i = 0; i < size(); i++) {
				index_type new_i = res.check_for_index(get_id(i), found);
				if (found) {
					res[new_i] = my_vector[i];
				}
			}
			swap(res);
		}
	}
	;

	/** Remap *this according to the index_to_variable_id_map iimap.
	 * Any variable not in iimap is disregarded.
	 *
	 * Returns true if all coefficients for variables
	 * not in dom are zero.
	 * Differs from reorder in that it doesn't throw. */
	virtual bool remap(const positional_vdomain& dom) {
		bool allzer = true;
		const index_to_variable_id_map_ptr& iimap = dom.get_index_to_variable_id_map();
		if (iimap != get_index_to_variable_id_map()) {
			vdom_vector res(iimap);
			bool found;
			for (unsigned int i = 0; i < size(); i++) {
				index_type new_i = res.check_for_index(get_id(i), found);
				if (found) {
					res[new_i] = my_vector[i];
				} else {
					if (allzer && !numeric::is_MEQ(my_vector[i], scalar_type(0)))
						allzer = false;
				}
			}
			swap(res);
		}
		return allzer;
	}
	;

	/** Remove the variables in vis, i.e., remove them from the domain. */
	virtual void remove_variables(const variable_id_set& vis) {
		index_to_variable_id_map_ptr iimap =
				get_index_to_variable_id_map()->get_map_with_ids_removed(vis);
		remap(iimap);
	}
	;

	/** Process the visitor v. */
	virtual void accept(vdom_vector_visitor<scalar_type>& v) const {
		v.visit(*this);
	}
	;

	/**
	 * Converts the elements of the vector to the type result_type.
	 * Uses the template function convert_element, which can be specialized to
	 * implement specific conversions between types.
	 */
	template<typename result_type> vdom_vector<result_type> convert_to() const {
		vector<result_type> v = my_vector.convert_to<result_type> ();
		vdom_vector<result_type> res(v, get_index_to_variable_id_map());
		return res;
	}
	;

	/** Prints the named vector.
	 * The default format is
	 * x1=a1,x2=a2, ... , xn=an.
	 */
	virtual void print(std::ostream& os) const;

	class output_format;
	static const output_format& get_output_format() {
		return my_output_format;
	}
	;
	static void set_output_format(output_format new_format) {
		my_output_format = new_format;
	}
	;

protected:
	friend vdom_vector<scalar_type>& operator+=<> (
			vdom_vector<scalar_type>& v1, const vdom_vector<scalar_type>& v2);
	friend vdom_vector<scalar_type>& operator*=<> (
			vdom_vector<scalar_type>& v, const scalar_type& c);
	friend vdom_vector<scalar_type>& operator/=<> (
			vdom_vector<scalar_type>& v, const scalar_type& c);
	/**
	 * Get coefficients a1,...,an as vector<scalar_type>.
	 */
	vector_type& get_impl_vector() {
		return my_vector;
	}
	;

	vector_type my_vector;
	static output_format my_output_format;
	static scalar_type my_zero; // to avoid constructing zero and being able to return a ref to const
	//static const scalar_type my_pos_one; // to avoid constructing zero and being able to return a ref to const
	//static const scalar_type my_neg_one; // to avoid constructing zero and being able to return a ref to const
};

/* need to declare static members */
template<typename scalar_type> typename vdom_vector<scalar_type>::output_format
		vdom_vector<scalar_type>::my_output_format;
template<typename scalar_type> scalar_type vdom_vector<scalar_type>::my_zero =
		scalar_type(0);

/*
 template <typename scalar_type> typename vdom_vector<scalar_type>::my_pos_one
 vdom_vector<scalar_type>::my_pos_one;
 template <typename scalar_type> typename vdom_vector<scalar_type>::my_neg_one
 vdom_vector<scalar_type>::my_neg_one;
 */

}

#include "math/vdom/vdom_vector_output.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/vdom/vdom_vector_utility.h"

#endif /*vdom_VECTOR_H_*/
