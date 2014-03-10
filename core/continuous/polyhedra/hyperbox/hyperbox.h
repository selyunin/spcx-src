/*
 * hyperbox.h
 *
 *  Created on: Nov 23, 2009
 *      Author: frehse
 */

// In order to break circular dependency with polyhedron_output.hpp,
// we apply a two-phase lock on the header
#if !defined(HYPERBOX_H_) || !defined(HYPERBOX_H__BODY)
//#ifndef HYPERBOX_H_
#define HYPERBOX_H_
#include "core/continuous/polyhedra/polyhedron.h"

#include <stdexcept>
#include <iostream>

#include "math/vdom/index_to_variable_id_map_provider.h"
#include "math/vdom/affine_map.h"
#include "math/vector.h"
#include "math/vector_operators.h"
#include "math/scalar_types/scalar_with_infinity.h"

#ifndef HYPERBOX_H__BODY
#define HYPERBOX_H__BODY

namespace continuous {

/** A hyperbox represented by its lower left and upper right corner.
 *
 * A hyperbox (l,u,iimap) defines the set
 *    l(x) <= x <= u(x)
 * for all x in iimap. The l(x) and u(x) are scalars extended with
 * infinity and NaN.
 *
 * In order to be consistent with polyhedra, a zero-dimensional
 * hyperbox can be empty (false) or not (true). */
template<typename scalar_type>
class hyperbox: public polyhedron<scalar_type> ,
		public index_to_variable_id_map_provider,
		public boost::enable_shared_from_this<hyperbox<scalar_type> > {
public:
	typedef hyperbox<scalar_type> my_type;
	typedef scalar_with_infinity<scalar_type> value_type;
	typedef math::vector<scalar_type> finite_point_type;
	typedef math::vector<value_type> point_type;
	typedef math::vdom_vector<scalar_type> vdom_vector_type;
	typedef typename point_type::size_type size_type;

	/** Create a universe hyperbox. */
	hyperbox();

	/** A hyperbox from point l to u, both defined over a variable domain.
	 *
	 * A point is a scalar_with_infinity, which can be a scalar,
	 * pos. or neg. infinity, or undefined. */
	hyperbox(point_type l, point_type u, const positional_vdomain& dom);

	/** A hyperbox from point l to u, both defined over iimap.
	 *
	 * A point is a scalar_with_infinity, which can be a scalar,
	 * pos. or neg. infinity, or undefined. */
	hyperbox(point_type l, point_type u, index_to_variable_id_map_ptr iimap);

	/** A hyperbox from finite points l to u, both defined over iimap.
	 *
	 * A finite point is a scalar number. */
	hyperbox(finite_point_type l, finite_point_type u,
			index_to_variable_id_map_ptr iimap);

	/** A hyperbox defined over iimap, with l and u undefined (NaN).
	 */
	explicit hyperbox(index_to_variable_id_map_ptr iimap);

	/** A hyperbox defined over domain, with l and u undefined (NaN).
	 */
	explicit hyperbox(const positional_vdomain& dom);

	/** Copy constructor */
	hyperbox(const hyperbox<scalar_type>& box);

	/** Assignment constructor */
	hyperbox<scalar_type>& operator=(const hyperbox<scalar_type>& box);

	virtual ~hyperbox();

	// ----------------------------------------
	// Continuous_set interface
	// ----------------------------------------

	/** Return a shared_ptr to *this. */
	virtual continuous_set::ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	virtual continuous_set::const_ptr get_const_ptr() const;

	virtual my_type* clone() const;

	virtual my_type* create_universe() const;
	virtual my_type* create_empty() const;

	virtual int get_memory() const;

	virtual unsigned int get_dim() const;

	virtual math::tribool is_empty() const;
	virtual math::tribool is_universe() const;

	virtual void embed_variables(const variable_id_set& id_set);
	virtual void
	existentially_quantify_variables(const variable_id_set& id_set);

	// Visitors are handled by polyhedron base class
	//	/** Accept a visitor. */
	//	virtual void
	//			accept(dispatching::dispatcher<continuous_set_typelist>& d) const;

	// ----------------------------------------
	// Support Function functions
	// ----------------------------------------

	/** Returns true if compute_support returns a support vector
	 * and false otherwise.
	 */
	virtual bool computes_support_vector() const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(
			const math::vdom_vector<scalar_type>& l,
			value_type& max_value,
			math::vdom_vector<value_type>& support_vec) const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(
			const math::vdom_vector<Rational>& l,
			Rational& max_value,
			math::vdom_vector<Rational>& support_vec, bool& is_empty,
			bool& is_bounded) const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<double>& l,
			double& max_value,
			math::vdom_vector<double>& support_vec, bool& is_empty,
			bool& is_bounded) const;

	// ----------------------------------------
	// Polyhedron functions
	// ----------------------------------------

	/** Get the constraints of *this. */
	virtual typename math::lin_constraint_system<scalar_type>::const_ptr
	get_constraints() const;

	/*! Adds the constraint \p c to \p *this.
	 */
	virtual void add_constraint(
			const math::lin_constraint<scalar_type> &c, bool check_redundancy = false);

	/** Remove all redundant constraints.
	 *
	 * Does nothing, since no redundant constraints.
	 * */
	virtual void remove_redundant_constraints();

	// ----------------------------------------
	// Hyperbox special functions
	// ----------------------------------------

	bool is_l_finite() const;

	bool is_u_finite() const;

	virtual bool is_finite() const;

	/** Compute the center *this.
	 *
	 * The calculation is done using scalar_with_infinity. If a variable is
	 * unbounded in one direction, the center is at infinity in that direction
	 * (positive or negative infinity).
	 * If a variable is infinity in both directions, the center is NOTANUMBER.
	 *
	 * Throws if *this is an empty box. */
	virtual point_type compute_center() const;

	/** Compute the center *this.
	 *
	 * Throws if *this is an empty box or one of the variables is unbounded. */
	virtual vdom_vector_type compute_finite_center() const;

	/** Returns the lower corner point. */
	virtual const point_type& get_l() const;

	/** Returns the lower corner point. */
	virtual const point_type& get_u() const;

	/** Returns the lower corner point.
	 *
	 * Throws if not finite.*/
	virtual vdom_vector_type get_finite_l() const;

	/** Returns the lower corner point.
	 *
	 * Throws if not finite.*/
	virtual vdom_vector_type get_finite_u() const;

	/** Reorder *this according to the domain dom.
	 *
	 * Throws if there is a nonzero coeffients for a variable
	 * that is not in the domain. If dom has new variables, the
	 * corresponding coefficients are set to zero.
	 * */
	virtual void reorder(const positional_vdomain& dom);

	/** Returns the Minkowski sum */
	hyperbox<scalar_type> operator+(const hyperbox<scalar_type>& h) const;

	/** Minkowski sum */
	hyperbox<scalar_type>& operator+=(const hyperbox<scalar_type>& h);

	/** Return an empty box over the variables in iimap, in canonic form. */
	static hyperbox<scalar_type> empty_box(index_to_variable_id_map_ptr iimap);

	/** Return an empty box over the variables in the domain, in canonic form. */
	static hyperbox<scalar_type> empty_box(const positional_vdomain& dom =
			positional_vdomain());


protected:
	/** All access to elements of l and u is via get and set methods,
	 * so that derived classes can use on-the-fly computation. */
	virtual const value_type& get_l(size_type i) const;

	virtual const value_type& get_u(size_type i) const;

	virtual void set_l(size_type i, value_type v);

	virtual void set_u(size_type i, value_type v);

	void set_zero_dim_empty();

private:
	/** The hyperbox is consistent if l[i] != +infty and u[i] != -infty.
	 */
	bool check_consistency();

	point_type my_l, my_u;
	bool my_zero_dim_empty;
};

}

#include "hyperbox.hpp"

#endif /* HYPERBOX_H__BODY */
#endif /* HYPERBOX_H_ */
