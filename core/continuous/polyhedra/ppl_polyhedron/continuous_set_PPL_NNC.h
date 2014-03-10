#ifndef CONTINUOUS_SET_PPL_NNC_H_
#define CONTINUOUS_SET_PPL_NNC_H_

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
#include "core/continuous/continuous_set_transforms/continuous_set_transform_declarations.h"
#include "core/continuous/polyhedra/polyhedron.h"
#include "math/vdom/index_to_variable_id_map_provider.h"
#include "core/continuous/polyhedra/ppl_polyhedron/extended_ppl.h"


namespace ppl_polyhedron {

/** A polyhedron implementation of continuous sets, based on the PPL NNC_Polyhedron class.
 *
 * @note Constraints are stored in the form a1 x1 + ... an xn + b [sign] 0,
 * where [sign] is ==, >= or >.
 */
class continuous_set_PPL_NNC : public continuous::polyhedron<Rational>,
	public index_to_variable_id_map_provider {
public:
	typedef boost::shared_ptr<continuous_set_PPL_NNC> ptr;
	typedef boost::shared_ptr<const continuous_set_PPL_NNC> const_ptr;
	//typedef enum {TEXTUAL, DOUBLE_CONSTRAINTS, DOUBLE_GENERATORS} output_format;
	typedef Parma_Polyhedra_Library::NNC_Polyhedron my_poly_type;
	//typedef continuous::polyhedron<Rational>::output_format output_format;

	// --------------------------------------------
	// Constructors
	// --------------------------------------------

	//! Constructor with the type of set as argument
	continuous_set_PPL_NNC();
	virtual ~continuous_set_PPL_NNC();

	/** Copy constructor. */
	//continuous_set_PPL_NNC(const continuous_set_PPL_NNC& p):_mypoly(p._mypoly), {};

	/**
	 Creates a unverse set of dimension equal to number of entries in index_to_id_map, containing the entire state space.
	 */
	continuous_set_PPL_NNC(const index_to_variable_id_map_ptr& pnew_map);

	/** Construct from a polyhedron and iimap. */
	continuous_set_PPL_NNC(const my_poly_type& poly, const index_to_variable_id_map_ptr& pnew_map);

	/** Creates an identical copy of *this. */
	virtual continuous_set_PPL_NNC* clone() const;

	/** Creates a set of dimension equal to number of entries in index_to_id_map, containing the entire state space.
	 * The resulting continuous_set_PPL_NNC uses the same index_to_variable_id_map as *this. */
	virtual continuous_set_PPL_NNC* create_universe() const;

	/** Creates an empty set of dimension equal to the number of entries in index_to_variable_id_map.
	 * The resulting continuous_set_PPL_NNC uses the same index_to_variable_id_map as *this.
	 */
	virtual continuous_set_PPL_NNC* create_empty() const;

	/** Creates a set of dimension equal to number of entries in index_to_id_map, containing the entire state space.
	 * The resulting continuous_set_PPL_NNC uses pnew_map as index_to_variable_id_map.
	 */
	virtual continuous_set_PPL_NNC* create_universe(
			const index_to_variable_id_map_ptr pnew_map) const;

	/** Creates an empty set of dimension equal to the number of entries in index_to_id_map.
	 * The resulting continuous_set_PPL_NNC uses pnew_map as index_to_variable_id_map.
	 */
	virtual continuous_set_PPL_NNC* create_empty(
			const index_to_variable_id_map_ptr pnew_map) const;

	// --------------------------------------------
	// Non-modifying non-semantic methods
	// --------------------------------------------

	/*! Returns the memory consumed by \p *this in bytes.
	 *  A utility function.
	 */
	virtual int get_memory() const;

	// --------------------------------------------
	// Non-modifying semantic methods
	// --------------------------------------------

	/*! Returns the dimension of the state space.
	 *  A utility function.
	 */
	virtual dimension_t get_dim() const;

	/*! \brief
	 Auxiliary set manipulation function.
	 Returns <CODE>true</CODE> if and only if \p *this is
	 empty.
	 */
	virtual math::tribool is_empty() const;

	/*! \brief
	 Auxiliary set manipulation function.
	 Returns <CODE>true</CODE> if and only if \p *this contains
	 the entire state space.
	 The default implementation constructs the universe, difference assigns \p *this and tests for emptiness.
	 */
	virtual math::tribool is_universe() const;

	/*! \brief
	 Auxiliary set manipulation function.
	 Returns <CODE>true</CODE> if and only if the intersection of \p *this and \p *ps is
	 empty.
	 The default implementation corresponds to the mathematical definition, but it may be overridden by a more efficient implementation.
	 */
	virtual math::tribool is_disjoint_from(const continuous::continuous_set_const_ptr& ps) const;

	/*! \brief
	 Auxiliary set manipulation function.
	 Returns <CODE>true</CODE> if and only if \p *this contains \p *ps.
	 The default implementation corresponds to difference assign followed by checking emtpiness, but it may be overridden by a more efficient implementation.
	 */
	virtual math::tribool contains(const continuous::continuous_set_const_ptr& ps) const;

	/** Returns true if compute_support returns a support vector
	 * and false otherwise.
	 */
	virtual bool computes_support_vector() const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value, math::vdom_vector<Rational>& support_vec, bool& is_empty,
			bool& is_bounded) const;
	virtual void compute_support(const math::vdom_vector<double>& l, double& max_value,
			math::vdom_vector<double>& support_vec, bool& is_empty, bool& is_bounded) const;

	// --------------------------------------------
	// Domain manipulation functions
	// --------------------------------------------

	/** Expand the state space to incorporate the variables in id_set in the sense of embedding.
	 If an id in id_set is already in the state space, perform the operation using id_set without id.
	 */
	virtual void embed_variables(const variable_id_set& id_set);

	/** Existential quantification over the variables in id_set.
	 If an id in id_set is not in the state space, perform the operation using id_set without id.
	 The implementation is responsible for updating _index_to_id_map of \p *this.
	 */
	virtual void existentially_quantify_variables(const variable_id_set& id_set);

	// --------------------------------------------
	// Modifying methods with 0 continuous_set argument
	// --------------------------------------------

	/*! \brief
	 A utility function.
	 Iterface to a (possibly time consuming) simplication operator.
	 The default implementation is to do nothing.
	 */
	virtual void simplify();

	// --------------------------------------------
	// Modifying methods with 1 continuous_set argument
	// --------------------------------------------

	/** Swaps \p *this and \p s. */
	virtual void swap(ptr s);

	/** Assigns to \p *this the intersection of \p *this and \p *p.
	 */
	virtual void intersection_assign(const continuous::continuous_set_const_ptr p);

	/** Assigns to \p *this the intersection of \p *this and \p *p.
	 */
	virtual void intersection_assign(const continuous_set_PPL_NNC& p);

	/** Assigns to \p *this the union of \p *this and \p *p.
	 Calls continuous_common function do_union_assign(this,p) to permit type conversion.
	 */
	virtual void union_assign(const continuous::continuous_set_const_ptr p);

	/** Assigns to \p *this the states of \p *this that are not in \p *p.
	 Calls continuous_common function do_difference_assign(this,p) to permit type conversion.
	 */
	virtual void difference_assign(const continuous::continuous_set_const_ptr ps);

	/** Assigns to \p *this an overapproximation of the states of \p *this that are not in \p *ps.
	 Used as a cheaper alternative to \p difference_assign.
	 */
	virtual void cheap_difference_assign(const continuous::continuous_set_const_ptr ps);

	/** Get the constraints of *this. */
	virtual math::lin_constraint_system<Rational>::const_ptr get_constraints() const;

	/*! Adds the constraint \p c to \p *this.
	 */
	virtual void add_constraint(const math::lin_constraint<Rational> &c, bool check_redundancy = false);

	/** Adds the constraint \p c to \p *this.
	 */
	virtual void add_constraint(const Parma_Polyhedra_Library::Constraint &c, bool check_redundancy = false);

	/** Removes all redundant constraints.

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

	/** Returns the constraints of the polyhedra as constant
	 * reference to Constraint_System
	 */
	virtual const Parma_Polyhedra_Library::Constraint_System & constraints() const;

	virtual void accept(dispatching::dispatcher<continuous::continuous_set_typelist>& d) const;

	/** Compute constant bound time elapse. */
	virtual void assign_transformation(
			const continuous::constant_bound_time_elapse_transform& t);
	/** Apply a reset function. */
	virtual void assign_transformation(
			const continuous::reset_function_transform& t);

protected:
	const Parma_Polyhedra_Library::NNC_Polyhedron& get_poly() const {
		return _mypoly;
	}
	;

//	/** Output as a stream of characters in textual format. */
//	virtual void print_textual(std::ostream& os) const;
//	/** Output as a stream of characters as a list of constraints with coefficients of type double. */
//	virtual void print_double_constraints(std::ostream& os) const;
//	/** Output as a stream of characters as a list of generators with coefficients of type double. */
//	virtual void print_double_generators(std::ostream& os) const;

	/** Re-index and embed variables such that *this contains all variables of \p p.
	 * Returns a map from indices in p to indices in *this. */
	index_to_index_bimap map_to_common_iimap(const index_to_variable_id_map_ptr& p);
	;

	/** Re-index and embed variables such that *this contains all variables of \p p.
	 * Also returns a pointer to a remapped version of \p p that has the same
	 * index_to_variable_id_map as *this. */
	ptr map_to_common_space(const const_ptr& p);
	;

	/** Re-index and embed variables such that *this contains all variables of \p p.
	 * Also returns a pointer to a remapped version of \p p that has the same
	 * index_to_variable_id_map as *this. If no remapping is necessary, then the
	 * returned PPL_NNC_const_ptr is equal to \p p, avoiding to clone. */
	const_ptr map_to_common_space_const(const continuous::continuous_set_const_ptr& p);
	;

//	/** Returns non-const pointers to remapped versions of p and q, cloning always both p and q.
//	 * If no remapping is necessary, mapped_p points to a clone of p and mapped_q to a clone of q. */
//	static void get_common_space_versions(const continuous::continuous_set_const_ptr p,
//			const continuous::continuous_set_const_ptr q, ptr& mapped_p, ptr& mapped_q);

	/** Returns const pointers to remapped versions of p and q, cloning only when necessary.
	 * If no remapping is necessary, mapped_p points to p and mapped_q to q. */
	static void get_common_space_versions_both_const(
			const continuous::continuous_set_const_ptr p,
			const continuous::continuous_set_const_ptr q, const_ptr& mapped_p,
			const_ptr& mapped_q);

private:
	Parma_Polyhedra_Library::NNC_Polyhedron _mypoly;
};

}

#endif /*CONTINUOUS_SET_PPL_NNC_H_*/
