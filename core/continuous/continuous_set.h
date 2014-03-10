#ifndef GUARD_continuous_set_h
#define GUARD_continuous_set_h
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

//#include <vector>
//#include <set>
//#include <map>
#include <iostream>
#include "boost/shared_ptr.hpp"
#include <boost/enable_shared_from_this.hpp>
#include "utility/shared_ptr_output.h"
#include "utility/printable.h"
#include "utility/tree_node.h"
#include "math/vdom/primed_variable_provider.h"
#include "core/continuous/continuous_set_declarations.h"
#include "utility/dispatching/double_dispatch.h"

namespace continuous {

class continuous_set_transform;
typedef boost::shared_ptr<const continuous_set_transform> continuous_set_transform_const_ptr;

typedef unsigned int dimension_t;

class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;

typedef tree::node continuous_set_predicate;

/**
 This class defines the base class for representations of continuous sets over a set of variables.

 Variables can be added (embedding the set in a state space of higher
 dimension) and removed (existential quantification).

 Intersection, union and difference operators might change the type
 of the object. E.g., the difference between two convex polyhedra might
 be a collection of convex polyhedra.
 Adding a new type should not interfere with old interfere with previous
 definitions. Therefore, all operators involving mixed types should
 be placed in the continuous_set_common class.

 A (continuous) relation is a continuous set over unprimed and singly-primed variables.
 We consider continuous sets as relations because any continuous set can be interpreted
 as a relation and vice versa. There does not seem to be a benefit from defining two
 seperate classes.
 */

typedef continuous_set relation;

class continuous_set: public virtual printable, public virtual primed_variable_provider,
public boost::enable_shared_from_this<continuous_set> {
public:
	typedef boost::shared_ptr<continuous_set> ptr;
	typedef boost::shared_ptr<const continuous_set> const_ptr;

	typedef dispatching::dispatcher<continuous_set_typelist> const_visitor;

	// --------------------------------------------
	/** \name Constructors
	 *  \{ */
	// --------------------------------------------

	/** Creates an identical copy of *this. */
	virtual continuous_set* clone() const = 0;

	/** Creates a set of dimension zero, containing the entire state space
	 * (equivalent to true). */
	virtual continuous_set* create_universe() const = 0;

	/** Creates an empty set of dimension zero (equivalent to false). */
	virtual continuous_set* create_empty() const = 0;

	/** Virtual Destructor. */
	virtual ~continuous_set() {
	}
	;

	/* \} */
	// --------------------------------------------
	/** \name Non-modifying non-semantic methods
	 *  \{ */
	// --------------------------------------------

	/** Return a shared_ptr to *this. */
	ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	const_ptr get_const_ptr() const;

	/** Returns the memory consumed by \p *this in bytes. */
	virtual int get_memory() const = 0;

	/* \} */
	// --------------------------------------------
	/** \name Non-modifying semantic methods
	 *  \{ */
	// --------------------------------------------

	/** Returns the dimension of the state space.
	 *  A utility function.
	 */
	virtual dimension_t get_dim() const = 0;

	/** Returns <CODE>true</CODE> if and only if \p *this is empty. */
	virtual math::tribool is_empty() const = 0;

	/** Returns <CODE>true</CODE> if and only if \p *this contains
	 the entire state space.
	 The default implementation constructs the universe, subtracts \p *this
	 and tests for emptiness.
	 */
	virtual math::tribool is_universe() const;

	/** Returns <CODE>true</CODE> if and only if the intersection of \p *this and \p *ps is
	 empty.
	 The default implementation corresponds to the mathematical definition, but it may be
	 overridden by a more efficient implementation.
	 */
	virtual math::tribool is_disjoint_from(const continuous_set_const_ptr& ps) const;

	/** Returns <CODE>true</CODE> if and only if \p *this contains \p *ps.
	 The default implementation corresponds to difference assign followed by
	 checking emptiness, but it may be overridden by a more efficient implementation.
	 */
	virtual math::tribool contains(const continuous_set_const_ptr& ps) const;

	/** Iterface to a (possibly time consuming) simplication operator.
	 The default implementation is to do nothing.
	 */
	virtual void simplify() {
	}
	;

	/* \} */
	// --------------------------------------------
	/** \name Domain manipulation
	 *  \{ */
	// --------------------------------------------

	/** Expand the state space to incorporate the variables in id_set in the sense of embedding.
	 If an id in id_set is already in the state space, perform the operation using id_set without id.
	 */
	virtual void embed_variables(const variable_id_set& id_set) = 0;

	/** \brief Existential quantification over the variables in id_set.
	 *
	 If an id in id_set is not in the state space, perform the operation using id_set without id.
	 */
	virtual void existentially_quantify_variables(const variable_id_set& id_set) = 0;

	/** \brief Project to the variables in id_set.
	 *
	 * Performs existential quantification over all variables except those in id_set.
	 */
	virtual void project_to_variables(const variable_id_set& id_set);

	/* \} */

	/** Returns the set in predicate form.*/
	virtual continuous_set_predicate::ptr get_predicate() const = 0;

	/**  Accept a const_visitor. A const_visitor must provide the function
	 * <code> void dispatch(const T* c) </code>
	 * for all derived classes of continuous_set that are listed in
	 * continuous_set_typelist (see continuous_set_declarations). */
	virtual void accept(const_visitor& d) const = 0;

};

}

#endif
