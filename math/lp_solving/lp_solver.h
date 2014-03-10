#ifndef lp_solver_H_
#define lp_solver_H_

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
#include "math/vdom/lin_constraint.h"
#include "math/vdom/lin_constraint_system.h"
//#include "support_function_provider.h"
#include <boost/shared_ptr.hpp>

namespace math {

/** An interface for linear programming solvers.
 * scalar_type must support the instantiations scalar_type(0) and scalar_type(1).
 */
template<typename scalar_type> class lp_solver {
public:
	typedef boost::shared_ptr<lp_solver<scalar_type> > ptr;
	typedef boost::shared_ptr<const lp_solver<scalar_type> > const_ptr;

	/** Virtual Destructor. */
	virtual ~lp_solver() {
	}
	;

	struct lp_result {
		scalar_type max_value;
		vdom_vector<scalar_type> support_vec;
		bool is_unsat;
		bool is_bounded;
		bool is_maximum;
	};

	/** Solve the linear program given by the objective function f under the constraints
	 * con. Formaly, compute max f(x) s.t. con(x).
	 * If f==0, check only satisfiability (support_vec,is_bounded,is_maximum undefined).
	 */
	virtual void maximize(const lin_expression<scalar_type>& f, const lin_constraint_system<
			scalar_type>& con, lp_result& res) = 0;

	/** Define the constraints *con for later (repeated) solving. */
	virtual void set_constraints(const lin_constraint_system<scalar_type>& con) = 0;

	/** Solve the linear program given by the objective function f under constraints
	 * defined previously with set_constraints(con). Formaly, compute max f(x) s.t. con(x).
	 * If f==0, check only satisfiability (support_vec,is_bounded,is_maximum undefined).
	 */
	virtual void maximize(const lin_expression<scalar_type>& f, lp_result& res) = 0;

	/** Get a string name for the solver. */
	virtual std::string get_solver_name() const = 0;
};

}

/** Stream output for lp_result. */
template<typename scalar_type> std::ostream& print(std::ostream& os,
const		typename math::lp_solver<scalar_type>::lp_result& res) {
			if (res.is_unsat)
			os << "unsatisfiable";
			else {
				if (res.is_bounded) {
					if (res.is_maximum)
					os << "max ";
					else
					os << "sup ";
					os << res.max_value << " at " << res.support_vec << "";
				} else {
					os << "unbounded";
				}
			}
			return os;
		}
;

#endif /*lp_solver_H_*/
