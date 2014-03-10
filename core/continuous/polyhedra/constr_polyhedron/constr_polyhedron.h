#ifndef constr_polyhedron_H_
#define constr_polyhedron_H_

/***************************************************************************
 *   Copyright (C) 2008 by Goran Frehse   *
 *   goran.frehse@imag.fr   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdexcept>
#include "core/continuous/polyhedra/polyhedron.h"
//#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "utility/calc_string.h"
#include "math/matrix.h"
#include "math/matrix_operators.h"
#include "math/lp_solving/lp_solver_user.h"

#include <iostream>

namespace continuous {

/** Forward declaration of classes used in header file. */
class continuous_set_transform;

typedef boost::shared_ptr<continuous_set_transform>
		continuous_set_transform_ptr;

math::tribool containment_test(const continuous_set_const_ptr &p1,
		const continuous_set_const_ptr &p2);

/**
 A polyhedron implementation of continuous sets, based on constraint
 representation with lp solver back-end. This class owns a set of constraints, on
 which all methods are based.
 */
template<typename scalar_type>
class constr_polyhedron: public polyhedron<scalar_type> ,
		public math::lp_solver_user<scalar_type> {
public:
	typedef constr_polyhedron<scalar_type> my_type;
	typedef math::lin_constraint_system<scalar_type> my_poly_type;
	typedef boost::shared_ptr<my_poly_type> my_poly_ptr;
	typedef boost::shared_ptr<my_type> ptr;
	typedef boost::shared_ptr<const my_type> const_ptr;
	typedef math::vdom_vector<scalar_type> point_type;

	/* --------------------------------------------
	 Constructors
	 -------------------------------------------- */

	/** Constructor, definition at the end of this file */
	constr_polyhedron();

	/** specify math::lp_solver */
	explicit constr_polyhedron(math::lp_solver<scalar_type> *s);

	/** Assign orig_poly to *this directly (without copying).  */
	explicit constr_polyhedron(my_poly_ptr orig_poly);

	/** Copy Constructor
	 *
	 * @remark Needs to call copy constructor of base classes, in particular of lp_solver_user. */
	constr_polyhedron(const my_type &orig_poly);

	/** Assignment
	 *
	 * @remark Needs to call assignment oeprator of base classes, in particular of lp_solver_user.*/
	constr_polyhedron<scalar_type>& operator=(const my_type &orig_poly);

	/** Destructor */
	virtual ~constr_polyhedron();

	/**
	 Creates a deep copy of *this.
	 */
	virtual my_type* clone() const;

	/*!
	 Creates a set of dimension zero, containing the entire state
	 space (true).
	 */
	virtual my_type* create_universe() const;

	/*!
	 Creates an empty set of dimension zero (false).
	 */
	virtual my_type* create_empty() const;

	/*!
	 Creates an empty set of dimension zero (false).
	 */
	static my_type empty_poly();


	/* --------------------------------------------
	 Non-modifying non-semantic methods
	 -------------------------------------------- */

	/** Returns the memory consumed by \p *this in bytes.
	 *  A utility function.
	 */
	virtual int get_memory() const;

	/* --------------------------------------------
	 Non-modifying semantic methods
	 -------------------------------------------- */

	/** Returns the variables in the constraints. */
	virtual const variable_id_set& get_variable_ids() const;

	/**
	 Returns the dimension of the state space. A utility function.
	 */
	virtual dimension_t get_dim() const;

	/**
	 Returns <CODE>true</CODE> if one of the constraints is unsatisfiable.
	 */
	virtual bool is_trivially_empty() const;

	/**
	 Returns <CODE>true</CODE> if and only if \p *this is empty.
	 */
	virtual math::tribool is_empty() const;

	/**
	 Tests if \p *this contains the entire state space.
	 This is true if and only if all constraints are satisfied by all
	 valuations.
	 */
	virtual math::tribool is_universe() const;

	/**
	 * Returns if a point is inside the polyhedron
	 */
	virtual math::tribool contains(const point_type& x) const;

	/**
	 * Pull the base class functions in, otherwise they get shadowed by the above definition.
	 */
	using polyhedron<scalar_type>::contains;

	/**
	 Removes some, but not necessarily all redundant constraints
	 With O(number of constraints) calls to the underlying solver and at most
	 one copying of underlying constraints to the solver
	 */
	virtual void simplify();

	/**
	 Removes all redundant constraints.

	 Since we handle equalities, redundancy can be hard to define.
	 We define redundancy like this:
	 Let A be the set of constraints, and c a constraint in A.
	 if A - {c} implies c,
	 then A := A - {c}.
	 else if (c is an equality e = b) and A - {c} implies c <= b
	 then A := A - {c} + e >= b
	 else if (c is an equality e = b) and A - {c} implies c >= b
	 then A := A - {c} + e <= b
	 */
	virtual void remove_redundant_constraints();

	/**
	 * Collapses inequalities to equality constraint when possible.
	 */
	virtual void collapse_inequalities();

	/** Returns true if compute_support returns a support vector
		 * and false otherwise.
		 */
	virtual bool computes_support_vector() const;

	/**
	 Computes the support function for a given vector l, i.e.,
	 the supremum of l*x over all x in *this.
	 If l has a variable that is not constrained in *this, then the result is
	 considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<Rational> &l,
			Rational & max_value, math::vdom_vector<Rational> &support_vec,
			bool & is_empty, bool & is_bounded) const;

	/**
	 Computes the support function for a given vector l, i.e.,
	 the supremum of l*x over all x in *this.
	 If l has a variable that is not constrained in *this, then the result is
	 considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<double>&l,
			double &max_value, math::vdom_vector<double>&support_vec,
			bool & is_empty, bool & is_bounded) const;

	virtual typename math::lin_constraint_system<scalar_type>::const_ptr
	get_constraints() const;

	/* --------------------------------------------
	 Methods for adding or removing ids
	 -------------------------------------------- */

	/**
	 Expand the state space to incorporate the variables in id_set
	 in the sense of embedding.
	 If an id in id_set is already in the state space, perform the
	 operation using id_set without id.
	 */
	virtual void embed_variables(const variable_id_set & id_set);

	/**
	 Existential quantification over the variables in id_set.
	 If an id in id_set is not in the state space, perform the
	 operation using id_set without id.
	 */
	virtual
	void existentially_quantify_variables(const variable_id_set &id_set);

	/* --------------------------------------------
	 Modifying methods with 1 continuous_set argument
	 -------------------------------------------- */

	/**
	 Swaps \p *this and \p s.
	 */
	virtual void swap(my_type& s);

	/**
	 Assigns to \p *this the intersection of \p *this and \p p.
	 */
	virtual void intersection_assign(const polyhedron<scalar_type> &p);

	/**
	 Adds the constraint \p c to \p *this.
	 */
	virtual void add_constraint(const math::lin_constraint<scalar_type> &c, bool check_redundancy = false);

	/** Returns the ids of the variables that are primed to degree \p prime_count. */
	virtual variable_id_set get_primed_variables(unsigned int prime_count) const;

	/** \name Methods changing the primedness of variables
	 *  \{ */

	/**
	 Set the primedness of the variables with primedness
	 of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int d, unsigned int p = 0);

	/**
	 Increase the primedness of the variables with primedness
	 of degree \p d by 1. If d is 0, increase all.
	 */
	virtual void increase_primedness(unsigned int d = 0);

	/**
	 Decrease the primedness of the variables with primedness
	 of degree \p d by 1.  If d is 0, decrease all.
	 */
	virtual void decrease_primedness(unsigned int d = 0);

	/* \} */

	my_poly_type & get_poly() const;

	/** Set the constraints to be a copy of new_poly. */
	void set_poly(const my_poly_type &new_poly);

	/** Set the constraints to be *new_poly. */
	void set_poly(my_poly_ptr new_poly);

	virtual void
	accept(dispatching::dispatcher<continuous_set_typelist> &d) const;

protected:
	/** Remove redundant constraints before vertice output. */
	virtual void print_double_generators(std::ostream& os) const;
	/** Remove redundant constraints before JVX output. */
	virtual void print_JVX(std::ostream& os) const;

private:
	math::tribool is_redundant(const math::lin_constraint<scalar_type> &con) const;

	virtual void maximize(const math::lin_expression<scalar_type> &f,
			typename math::lp_solver<scalar_type>::lp_result &res) const;

	my_poly_ptr my_poly;
	/* dirty bit for updating solver w/ constraints */
	mutable bool my_up_to_date;

}; /* class constr_polyhedron */

} // namespace continuous

/** For calc_string as scalar_type, the computations (for enumerating vertices etc.) don't
 * work, so we redirect to a different routine. */
std::ostream& operator<<(std::ostream &os, const continuous::constr_polyhedron<
		calc_string>::ptr &p);

#include "constr_polyhedron.hpp"
//#include "constr_polyhedron_operators.h"

#endif /*constr_polyhedron_H_ */
