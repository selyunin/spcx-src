#ifndef lin_constraint_H_
#define lin_constraint_H_

//#include "global/global_types.h" // for __float128 operator<<

#include <boost/shared_ptr.hpp>
#include "utility/operator_enums.h"
#include "utility/shared_ptr_user.h"
#include "math/vdom/lin_expression.h"
#include "math/numeric/approx_comparator.h"
#include "math/numeric/comp.h"

namespace math {

template<typename scalar_type> class lin_constraint;

/** Visitor for processing lin_constraint objects */
template<typename scalar_type> class lin_constraint_visitor: public lin_expression_visitor<
		scalar_type> {
public:
	virtual ~lin_constraint_visitor() {
	}
	;
	virtual void visit(const lin_constraint<scalar_type>& c) {
		lin_constraint_prologue(c);
		c.get_l().accept(*this);
		visit_sign(c, c.get_sign());
		lin_constraint_epilogue(c);
	}
	;
	/** The prologue is called before processing the linear expression. */
	virtual void lin_constraint_prologue(const lin_constraint<scalar_type>& c) {
	}
	;
	/** The visit_sign is called for the sign */
	virtual void visit_sign(const lin_constraint<scalar_type>& c,
			const comparison_operator& s) {
	}
	;
	/** The epilogue is called after processing the linear expression. */
	virtual void lin_constraint_epilogue(const lin_constraint<scalar_type>& c) {
	}
	;
};

/** Represents a linear constraint of the form
 * a1 x1 + a2 x2 + ... + an xn + b [sign] 0,
 * where [sign] is in {<,<=,>,>=,==}.
 * The canonic form is
 *   a1 x1 + a2 x2 + ... + an xn + b <= 0.
 * Extends the class lin_expression by adding a sign.
 * */
template<typename scalar_type> class lin_constraint: public primed_variable_provider,
		public printable,
		public ptr_interface<lin_constraint<scalar_type> > {
public:
	typedef boost::shared_ptr<lin_constraint<scalar_type> > ptr;
	typedef boost::shared_ptr<const lin_constraint<scalar_type> > const_ptr;
	typedef comparison_operator sign;

	/* The default constraint is 0<=0 (satisfiable). */
	lin_constraint() :
		my_l(), my_s(LE) {
	}
	;
	/** Construct a constraint a^Tx+b s 0, where s denotes the sign. */
	lin_constraint(const vdom_vector<scalar_type>& v, const scalar_type& b, sign s) :
		my_l(lin_expression<scalar_type>(v,b)), my_s(s) {
	}
	;
	lin_constraint(const lin_expression<scalar_type>& l, sign s) :
		my_l(l), my_s(s) {
	}
	;
	virtual ~lin_constraint() {
	}
	;

	lin_expression<scalar_type>& get_l() {
		return my_l;
	}
	;

	const lin_expression<scalar_type>& get_l() const {
		return my_l;
	}
	;

	/** Returns the normal vector of the constraints.
	 *
	 * The normal vector points toward the outside of the
	 * halfspace (for inequalities).
	 */
	vdom_vector<scalar_type> get_normal() const {
		return get_canonic_l().get_vdom_vec();
	}
	;

	/** Get the coefficient of variable id to the value x.
	 * If id is not in the domain of *this, return zero.
	 */
	const scalar_type& get_coeff_with_id(const variable_id& id) const {
		return my_l.get_coeff_with_id(id);
	}
	;
	/**
	 * Get inhomogeneous coefficient b.
	 */
	const scalar_type& get_inh_coeff() const {
		return my_l.get_inh_coeff();
	}
	;
	/** Returns true if the constraint is in
	 * canonic form.
	 */
	bool is_canonic() const {
		return !is_GT_or_GE(get_sign());
	}
	/** Get the coefficients brought to canonic form,
	 * i.e., multiplied with -1 if necessary.
	 */
	lin_expression<scalar_type> get_canonic_l() const {
		if (!is_canonic())
			return -get_l();
		else
			return get_l();
	}
	;
	/** Get the coefficient of variable id to the value x brought to canonic form,
	 * i.e., multiplied with -1 if necessary.
	 * If id is not in the domain of *this, return zero.
	 */
	scalar_type get_canonic_coeff_with_id(const variable_id& id) const {
		if (!is_canonic())
			return -get_coeff_with_id(id);
		else
			return get_coeff_with_id(id);
	}
	;
	/** Get the inh. coefficient b brought to canonic form,
	 * i.e., multiplied with -1 if necessary.
	 */
	scalar_type get_canonic_inh_coeff() const {
		if (!is_canonic())
			return -get_inh_coeff();
		else
			return get_inh_coeff();
	}
	;
	lin_constraint<scalar_type> get_canonical() const {
		if(!is_canonic())
			return lin_constraint<scalar_type>(get_canonic_l(), get_canonic_sign());
		else
			return lin_constraint<scalar_type>(my_l, my_s);
	}
	const sign& get_sign() const {
		return my_s;
	}
	;
	sign get_canonic_sign() const {
		if (!is_canonic())
			return get_flipped_sign(my_s);

		else
			return my_s;
	}
	;
	/** Return true if the sign is ==, and false otherwise. */
	bool is_equality() const {
		return is_EQ(my_s);
	}
	/** Return false if the sign is ==, and true otherwise. */
	bool is_inequality() const {
		return !is_equality();
	}
	/** Return true if the sign is < or >, and false otherwise. */
	bool is_strict_inequality() const {
		return is_strict(my_s);
	}
	/** Return true if the constraint is satisfiable, i.e., there is a
	 * valuation that satisfies the constraint.
	 * This is the case if and only if the homogeneous coefficients are
	 * not all zero or the evaluation of the constraint based on the inhomogeneous
	 * coefficient returns true. */
	math::tribool is_satisfiable() const {
		if (get_l().is_homogeneous_coeffs_zero())
			return evaluate_inhomogeneous_coefficient();
		else
			return true; // it is satisfiable
	}
	;
	/** Return true if the constraint is always satisfied, i.e., there any
	 * valuation satisfies the constraint.
	 * This is the case if and only if the homogeneous coefficients are
	 * all zero and the evaluation of the constraint based on the inhomogeneous
	 * coefficient returns true. */
	math::tribool is_always_satisfied() const {
		if (get_l().is_homogeneous_coeffs_zero())
			return evaluate_inhomogeneous_coefficient();
		else
			return false; // it is not always satisfied
	}
	;
	/** Evaluate the constraint based on the inhomogeneous coeffiecient, i.e.,
	 * consider the homogeneous coefficients to be zero.
	 *
	 *
	 * Uses numeric approx comparisons. */
	math::tribool evaluate_inhomogeneous_coefficient() const {
//std::cout << *this << " eval inh: " << tribool_string(eval_with_sign(get_inh_coeff(),scalar_type(0))) << std::endl;
		return eval_with_sign(get_inh_coeff(),scalar_type(0));
	}
	;
	/** Returns an unsatisfiable constraint, 1<=0. */
	static lin_constraint<scalar_type> unsatisfiable_constraint() {
		lin_expression<scalar_type> l;
		l.set_inh_coeff(scalar_type(1));
		return lin_constraint<scalar_type> (l, LE);
	}
	;

	/** Returns a*v.
	 */
	scalar_type normal_eval(const vdom_vector<scalar_type>& v) const {
		return scalar_product(get_normal(),v);
	}
	;

	/** Returns a*v + b.
	 */
	scalar_type eval(const vdom_vector<scalar_type>& v) const {
		return scalar_product(get_normal(),v) + get_canonic_inh_coeff();
	}
	;

	/** Returns true if the constraint is active on the point defined
	 * by v, i.e., whether a*v+b==0.
	 *
	 * For approx. scalar types, this computes whether a*v and -b are
	 * approx. equal.*/
	bool is_active(const vdom_vector<scalar_type>& v) const {
		const vdom_vector<scalar_type>& u = get_l().get_vdom_vec();
//std::cout << *this << " is active on " << v << ": if " << scalar_product(u,v) << " equals " << -get_inh_coeff() << std::endl;
		return math::numeric::is_MEQ(scalar_product(u,v),-get_inh_coeff());
	}
	;


	/** Returns true if the constraint is satisfied on the point defined
	 * by v, i.e., whether a*v+b sign 0.
	 *
	 * For approx. scalar types, this computes whether a*v and -b are
	 * approx. smaller, equal or larger (according to the sign).
	 * The result is conservative in the sense of "maybe" satisfied. */
	math::tribool is_satisfied(const vdom_vector<scalar_type>& v) const {
		const vdom_vector<scalar_type>& u = get_l().get_vdom_vec();
		scalar_type x=scalar_product(u,v);
		scalar_type y=-get_inh_coeff();
		// the constraint is satisfied if x < y, x <= y, respectively x == y
//std::cout << "is_active:" << scalar_product(u,v) << " vs " << -get_inh_coeff() << std::endl;
		return eval_with_sign(x,y);
	}
	;

	/** Zero all coefficients.
	 *  */
	void clear() {
		my_l.clear();
	}
	;

	/** Change the sign to new_sign. */
	void set_sign(const sign& new_sign) {
		my_s = new_sign;
	}
	;
	/** Set the coefficient of variable id to the value x.
	 * Note that this might be slow due to internal reassignment of the iimap.
	 */
	void set_coeff_with_id(const variable_id& id, const scalar_type& s) {
		my_l.set_coeff_with_id(id, s);
	}
	;
	/** Set inhomogeneous coefficient b.
	 */
	void set_inh_coeff(const scalar_type& b) {
		my_l.set_inh_coeff(b);
	}
	;

	/** Make the constraint closed. */
	void closure_assign() {
		my_s = get_closed_sign(my_s);
	}
	;

	/** Prints the constraint in the format
	 * a1 x1 + a2 x2 + ... + an xn + b <sign> 0
	 * Nonzero elements are omitted.
	 */
	virtual void print(std::ostream& os) const {
		/* print the linear constraint */

		// if the first nonzero coeff is negative, reverse the sign
		bool flip_sign = false;
		if (my_l.size() > 0) {
			unsigned int i = 0;
			while (i < my_l.size() && my_l[i] == scalar_type(0)) {
				++i;
			}
			if (i < my_l.size())
				flip_sign = my_l[i] < scalar_type(0);
		}
		if (flip_sign) {
			lin_constraint c(-my_l, get_flipped_sign(my_s));
			c.print(os);
		} else {
			my_l.print_hom(os);
			os << " " << sign_string(my_s) << " ";
			scalar_type c = -my_l.get_inh_coeff();
			os << c;
		}
	}

	/** Returns the ids of all variables over which the set is defined. */
	virtual const variable_id_set& get_variable_ids() const {
		return my_l.get_variable_ids();
	}
	;

	/** Returns the ids of the variables that are primed to degree \p prime_count. */
	virtual variable_id_set get_primed_variables(unsigned int prime_count) const {
		return my_l.get_primed_variables(prime_count);
	}
	;

	/** Remove the variables in vis, i.e., remove them from the domain. */
	virtual void remove_variables(const variable_id_set& vis) {
		my_l.remove_variables(vis);
	}
	;

	/**
	 * Converts the coefficients of the constraint to the type result_type.
	 * Uses the template function convert_element, which can be specialized to
	 * implement specific conversions between types.
	 */
	template<typename result_type> lin_constraint<result_type> convert_to() const {
		lin_constraint<result_type> res(
				my_l.template convert_to<result_type> (), my_s);
		return res;
	}
	;

	/** Process the visitor v. */
	virtual void accept(lin_constraint_visitor<scalar_type>& v) const {
		v.visit(*this);
	}
	;

	/** Returns a constraint without variables that is unsatisfiable. */
	static lin_constraint<scalar_type> zero_dim_false() {
		lin_constraint<scalar_type> l = zero_dim_true();
		l.set_inh_coeff(scalar_type(1));
		return l;
	}
	;

	/** Returns a constraint without variables that is satisfied
	 * for any variable valuation. */
	static lin_constraint<scalar_type> zero_dim_true() {
		return lin_constraint<scalar_type> ();
	}
	;

	// --------------------------------------------
	/** \name Methods changing the primedness of variables
	 *  \{ */
	// --------------------------------------------

	/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int d, unsigned int p = 0) {
		my_l.reassign_primedness(d, p);
	}
	;

	/** Increase the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int d = 0) {
		my_l.increase_primedness(d);
	}
	;

	/** Decrease the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int d = 0) {
		my_l.decrease_primedness(d);
	}
	;

	/* \} */
	// --------------------------------------------
protected:
	math::tribool eval_with_sign(const scalar_type& x, const scalar_type& y) const {
		using namespace math::numeric;
		switch (my_s) {
		case LT:
			return is_LT(x, y);
		case LE:
			return is_LE(x, y);
		case EQ:
			return is_EQ(x, y);
		case GE:
			return is_GE(x, y);
		case GT:
			return is_GT(x, y);
		default:
			throw std::runtime_error("eval_with_sign: unknown sign");
			return indeterminate();
		}
	}



	lin_expression<scalar_type> my_l;
	sign my_s;
};

}

#include "math/vdom/lin_constraint_operators.h"

#endif /*lin_constraint_H_*/
