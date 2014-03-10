/*
 * finite_finite_hyperbox.h
 *
 *  Created on: Jan 17, 2011
 *      Author: frehse
 */

// In order to break circular dependency with polyhedron_output.hpp,
// we apply a two-phase lock on the header
#if !defined(FINITE_HYPERBOX_H_) || !defined(FINITE_HYPERBOX_H__BODY)
//#ifndef HYPERBOX_H_
#define FINITE_HYPERBOX_H_
#include "core/continuous/polyhedra/polyhedron.h"

#include <stdexcept>
#include <iostream>

#include "math/vdom/index_to_variable_id_map_provider.h"
#include "math/vdom/affine_map.h"
#include "math/vector.h"
#include "math/vector_operators.h"

#ifndef FINITE_HYPERBOX_H__BODY
#define FINITE_HYPERBOX_H__BODY

namespace continuous {

/** A finite_hyperbox represented by its center and a generator.
 *
 * A finite_hyperbox (c,g,dom) defines the set x given component-wise as
 *    x_i = c_i + u_i*g_i, -1 <= u_i <= 1
 * for all x_i in dom.
 *
 * A finite_hyperbox can not be universe, nor zero-dimensional.
 * An exception is throws if such objects are requested
 * through calls of the base class interface.
 *
 * A finite_hyperbox is empty iff there is a g_i<0. */
template<typename scalar_type>
class finite_hyperbox: public polyhedron<scalar_type> ,
		public index_to_variable_id_map_provider,
		public boost::enable_shared_from_this<finite_hyperbox<scalar_type> > {
public:
	typedef boost::shared_ptr<finite_hyperbox<scalar_type> > ptr;
	typedef boost::shared_ptr<const finite_hyperbox<scalar_type> > const_ptr;

	typedef finite_hyperbox<scalar_type> my_type;
	typedef scalar_type value_type;
	typedef math::vector<scalar_type> point_type;
	typedef math::vdom_vector<scalar_type> vdom_vector_type;
	typedef typename point_type::size_type size_type;

	/** Inconsistent constructor */
	finite_hyperbox();

	/** A finite_hyperbox with center c and generator g, both defined over dom.
	 */
	finite_hyperbox(point_type c, point_type g, positional_vdomain dom);

	/** A finite_hyperbox with center c and generator g, both defined over the same domain.
	 */
	finite_hyperbox(vdom_vector_type c, vdom_vector_type g);

	/** Copy constructor */
	finite_hyperbox(const finite_hyperbox<scalar_type>& box);

	/** Assignment constructor */
	finite_hyperbox<scalar_type>& operator=(
			const finite_hyperbox<scalar_type>& box);

	// ----------------------------------------
	// Continuous_set interface
	// ----------------------------------------

	/** Return a shared_ptr to *this. */
	continuous_set::ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	continuous_set::const_ptr get_const_ptr() const;

	my_type* clone() const;

	my_type* create_universe() const;
	my_type* create_empty() const;

	int get_memory() const;

	unsigned int get_dim() const;

	math::tribool is_empty() const;
	math::tribool is_universe() const;

	void embed_variables(const variable_id_set& id_set);
	void
	existentially_quantify_variables(const variable_id_set& id_set);

	// Visitors are handled by polyhedron base class
	//	/** Accept a visitor. */
	//	 void
	//			accept(dispatching::dispatcher<continuous_set_typelist>& d) const;

	// ----------------------------------------
	// Support Function functions
	// ----------------------------------------

	/** Returns true if compute_support returns a support vector
	 * and false otherwise.
	 */
	bool computes_support_vector() const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	void
	compute_support_impl(const math::vdom_vector<scalar_type>& l,
			value_type& max_value, math::vdom_vector<value_type>& support_vec,
			bool& is_empty, bool& is_bounded) const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value, math::vdom_vector<Rational>& support_vec,
			bool& is_empty, bool& is_bounded) const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	void compute_support(const math::vdom_vector<double>& l, double& max_value,
			math::vdom_vector<double>& support_vec, bool& is_empty,
			bool& is_bounded) const;

	// ----------------------------------------
	// Polyhedron functions
	// ----------------------------------------

	/** Get the constraints of *this. */
	typename math::lin_constraint_system<scalar_type>::const_ptr
	get_constraints() const;

	/*! Adds the constraint \p c to \p *this.
	 *
	 * Throws if the constraint is not a box constraint.
	 */
	void add_constraint(const math::lin_constraint<scalar_type> &c,
			bool check_redundancy = false);

	/** Remove all redundant constraints.
	 *
	 * Does nothing, since no redundant constraints.
	 * */
	virtual void remove_redundant_constraints();


	// ----------------------------------------
	// Hyperbox special functions
	// ----------------------------------------

	/** Returns the center. */
	const point_type& get_c() const;

	/** Returns the generator. */
	const point_type& get_g() const;

	/** Returns the center with domain.
	 */
	vdom_vector_type get_c_dom() const;

	/** Returns the generator with domain.
	 */
	vdom_vector_type get_g_dom() const;

	/** Returns the lower bound */
	point_type lower() const;

	/** Returns the upper bound */
	point_type upper() const;

	/** Check if point is contained in the hyperbox
	 *
	 * @note Variables not in the domain of the hypebox are
	 * disregarded. I.e., the hyperbox is embedded in the
	 * common space in the sense of extending to infinity. */
	bool contains(const vdom_vector_type& p) const;

	/** Reorder *this according to the domain dom.
	 *
	 * Throws if there is a nonzero coeffients for a variable
	 * that is not in the domain. If dom has new variables, the
	 * corresponding coefficients are set to zero.
	 * */
	void reorder(const positional_vdomain& dom);

	/** Comparison */
	bool operator==(const finite_hyperbox<scalar_type>& h) const;

	/** Minkowski sum */
	finite_hyperbox<scalar_type> operator+(const finite_hyperbox<scalar_type>& h) const;

	/** Minkowski sum */
	finite_hyperbox<scalar_type>& operator+=(const finite_hyperbox<scalar_type>& h);

	/** Make the box empty. */
	void set_empty();

	/** Return an empty box over the variables in iimap, in canonic form. */
	static finite_hyperbox<scalar_type> empty_box(const positional_vdomain& dom);

private:
	/** Canonicalize by making g[0] neg if empty.
	 */
	void canonicalize();

	point_type my_c, my_g;
};

/** Specialization to avoid unnecessary conversions */
template<>
inline void finite_hyperbox<double>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const;

}

#include "finite_hyperbox.hpp"
#include "finite_hyperbox_operators.h"

#endif /* FINITE_HYPERBOX_H__BODY */
#endif /* FINITE_HYPERBOX_H_ */
