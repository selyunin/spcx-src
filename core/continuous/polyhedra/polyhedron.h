#ifndef polyhedron_H_
#define polyhedron_H_

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
#include <boost/shared_ptr.hpp>
#include "utility/printable.h"
#include "math/vdom/lin_constraint.h"
#include "math/vdom/lin_constraint_system.h"
#include "core/continuous/support_function_provider.h"
//#include "continuous_set.h"

namespace continuous {

/** An interface for polyhedron implementations of continuous sets.
 * The polyhedron is given by a collection of linear constraints.
 * scalar_type must support the instantiations scalar_type(0) and scalar_type(1).
 */
template<typename scalar_type> class polyhedron: public support_function_provider
//public virtual printable
//public ptr_interface<polyhedron<scalar_type> >
{
public:
	typedef enum {
		TEXTUAL, DOUBLE_CONSTRAINTS, DOUBLE_GENERATORS, JVX
	} output_format;

	typedef boost::shared_ptr<polyhedron<scalar_type> > ptr;
	typedef boost::shared_ptr<const polyhedron<scalar_type> > const_ptr;

	// --------------------------------------------
	// Constructors
	// --------------------------------------------

	/** Constructor */
	polyhedron() {
	}
	;

	/** Virtual Destructor.
	 */
	virtual ~polyhedron() {
	}
	;

	/** Virtual copy constructor.	 */
	virtual polyhedron<scalar_type>* clone() const = 0;

	/** Creates a set of dimension zero, containing the entire state space
	 * (equivalent to true). */
	virtual polyhedron<scalar_type>* create_universe() const = 0;

	/** Creates an empty set of dimension zero (equivalent to false). */
	virtual polyhedron<scalar_type>* create_empty() const = 0;

	/** Get the constraints of *this. */
	virtual typename math::lin_constraint_system<scalar_type>::const_ptr
	get_constraints() const = 0;

	/*! Adds the constraint \p c to \p *this.
	 *
	 * If check_redundancy is true, checks whether c is redundant,
	 * and only adds it if it is not.
	 * @attention Adding more than one constraints with redundancy check
	 * does not imply that there are no redundant constraints.
	 * A previously added constraint may be made redundant by subsequent
	 * added constraint. To ensure non-redundancy, use the function remove_redundant_constraints().
	 */
	virtual void add_constraint(
			const math::lin_constraint<scalar_type> &c, bool check_redundancy = false) = 0;

	/** Add the constraints in the set con_set.
	 *
	 * If check_redundancy is true, checks whether c is redundant,
	 * and only adds it if it is not.
	 * @attention Adding more than one constraints with redundancy check
	 * does not imply that there are no redundant constraints.
	 * A previously added constraint may be made redundant by subsequent
	 * added constraint. To ensure non-redundancy, use the function remove_redundant_constraints().
	 */
	virtual void add_constraints(const math::lin_constraint_system<
			scalar_type>& con_set, bool check_redundancy = false) {
		for (typename math::lin_constraint_system<scalar_type>::const_iterator
				i = con_set.begin(); i != con_set.end(); ++i) {
			add_constraint(*i, check_redundancy);
		}
	}
	;

	/** Add the constraints in the set con_set.
	 */
	virtual void add_constraints(
			const typename math::lin_constraint_system<scalar_type>::const_ptr& con_set) {
		add_constraints(*con_set);
	}
	;

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
	virtual void remove_redundant_constraints() = 0;

	/** Returns the set in predicate form.
	 *
	 * The default implementation obtains the constraints and constructs the
	 * predicates from them.*/
	virtual continuous_set_predicate::ptr get_predicate() const {
		throw std::runtime_error("missing implementation get_predicate");
		return continuous_set_predicate::ptr();
	}
	;

	virtual void accept(dispatching::dispatcher<continuous_set_typelist>& d) const {
		d.dispatch(this);
	}
	;

	/** Set the output format according to \p s:
	 * - TEXTUAL : in textual format (readable by the input parser)
	 * - DOUBLE_CONSTRAINTS : output the constraints as a sequence of double values separated by spaces, one constraint per line.
	 *   I.e., the constraint \f$a_0 x_0 + a_1 x_1 + ... + a_{n-1} x_{n-1} + b \geq 0\f$ is output as
	 *   "a0 a1 a2 ... an b".
	 *   There is no distinction between strict and nonstrict constraints. Equalities are converted to two nonstrict inequalities.
	 *   If the
	 *    polyhedron is universal, the list is empty.
	 * - DOUBLE_GENERATORS : output the vertices as a sequence of double values separated by spaces, one generator per line.
	 *   I.e., the vertice \f$(x_0,x_1, ... ,x_{n-1})\f$ is output as
	 *   "x0 x1 ... xn".
	 *   If the polyhedron is 2-dimensional (a polygon), the vertices are ordered counterclockwise, and the polygon is closed by
	 *   repeating the first point at the end. This is necessary for output with the tool "graph".
	 *   Rays are ignored in this format, i.e., only bounded polyhedra can be printed correctly.
	 *   Note that due to the internal representation of the PPL vertices can appear twice (as closure points).
	 *
	 * The default format is TEXTUAL.
	 */
	static void set_output_format(output_format f) {
		my_output_format = f;
	}
	;

	/** Get the current print format. Useful, e.g., for saving the current format to restore it later. */
	static output_format get_output_format() {
		return my_output_format;
	}
	;

	/**
	 * Output as a stream of characters. The format can be set using \p set_print_format.
	 */
	virtual void print(std::ostream& os) const {
		if (get_output_format() == TEXTUAL) {
			print_textual(os);
		} else if (get_output_format() == DOUBLE_CONSTRAINTS) {
			print_double_constraints(os);
		} else if (get_output_format() == DOUBLE_GENERATORS) {
			print_double_generators(os);
		} else if (get_output_format() == JVX) {
			print_JVX(os);
		} else {
			throw std::runtime_error("unknown polyhedron output format");
		}
	}
	;

	virtual void print_textual(std::ostream& os) const;
	virtual void print_double_constraints(std::ostream& os) const;
	virtual void print_double_generators(std::ostream& os) const;
	virtual void print_JVX(std::ostream& os,variable_id_list vil=variable_id_list()) const;

private:
	static output_format my_output_format;
};

template<typename s> typename polyhedron<s>::output_format
		polyhedron<s>::my_output_format;

}

//#include "core/continuous/polyhedra/polyhedron_operators.h"
#include "core/continuous/polyhedra/polyhedron_output.h"

#endif /*polyhedron_H_*/
