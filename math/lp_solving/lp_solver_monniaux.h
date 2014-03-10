#ifndef LP_SOLVER_MONNIAUX_H_
#define LP_SOLVER_MONNIAUX_H_

#define BOOST_VECTORS

inline double numericConversion(const Rational& q) {
	return q.get_double();
}

#include "monniaux/fastSimplex.hh"
#include "math/lp_solving/lp_solver.h"
#include <boost/shared_ptr.hpp>

namespace math {

/** A linear programming solver based on the simplex algorithm of David Monniaux.
 * scalar_type must support the instantiations scalar_type(0) and scalar_type(1).
 */
template<typename scalar_type> class lp_solver_monniaux: public lp_solver<scalar_type> {
public:
	typedef boost::shared_ptr<lp_solver_monniaux<scalar_type> > ptr;
	typedef boost::shared_ptr<const lp_solver_monniaux<scalar_type> > const_ptr;

	/** Virtual Destructor. */
	virtual ~lp_solver_monniaux() {
	}
	;

	/* Convert a lin_constraint to LinearForm using the indices in iimap. */
	LinearForm<scalar_type> convert_to_LinearForm(const lin_constraint<scalar_type>& c,
			const std::vector<std::string>& variableNames,
			index_to_variable_id_map::const_ptr iimap) {
		/* To do : conversion according to the sign of the constraint! */
		LinearForm<scalar_type> lf(variableNames);
		for (unsigned int i = 0; i < c.get_l().size(); ++i) {
			index_type ind = iimap->get_index(c.get_l().get_iimap()->get_id(i));
			lf[ind] = c.get_l()[i];
		}
		return lf;
	}
	;

	/** Solve the linear program given by the objective function f under the constraints
	 * con. Formaly, compute max f(x) s.t. con(x).
	 */
	virtual void maximize(const math::lin_expression<scalar_type>& f, //
const			typename lin_constraint_system<scalar_type>::const_ptr& con,
			typename lp_solver<scalar_type>::lp_result& res) {
				/* get all variables in the constraints */
				variable_id_set actual_variables;
				for (typename lin_constraint_system<scalar_type>::const_iterator it=
						con->begin(); it !=con->end(); ++it) {
					variable_id_set vis=(*it).get_variable_ids();
					actual_variables.insert(vis.begin(), vis.end());
				}

				/* get the total number of variables
				 *   = actual variables + number of inequality constraints. */

				/* count inequality constraints */
				unsigned int ineq_count=0;
				for (typename lin_constraint_system<scalar_type>::const_iterator it=
						con->begin(); it !=con->end(); ++it) {
					if ((*it).get_sign()!=lin_constraint<scalar_type>::EQ)
					++ineq_count;
				}

				unsigned int actual_variable_count = actual_variables.size();
				unsigned int variable_count = actual_variable_count+ineq_count;

				/* Get a vector of variable names
				 */
				std::vector<std::string> variableNames(variable_count);
				{
					unsigned int i=0;
					for (variable_id_set::const_iterator it=actual_variables.begin(); it
							!= actual_variables.end(); ++it) {
						variableNames[i]=variable::get_name(*it);
						++i;
					}
					/* fill the rest */
					for (; i<variable_count; ++i) {
						variableNames[i]="aux"+int2string(i);
					}
				}

				/* Instantiate solver. */
				SimplexF active_solver(variableNames, actual_variable_count);

				/* Add constraints. */
				/* Prepare mapping the variables to common indices. */
				index_to_variable_id_map::const_ptr iimap =
				index_to_variable_id_map::get_map_with_ids(actual_variables);

				for (typename lin_constraint_system<scalar_type>::const_iterator it=
						con->begin(); it !=con->end(); ++it) {
					LinearForm<scalar_type> lf=convert_to_LinearForm(*it, variableNames,
							iimap);

					/* To do : How to get i? */
					unsigned int i;
					active_solver.assignBasic(i, lf);
				}

				/* Assign const function */

				/* Solve optimization problem */

				/* Define result */

			}
			;

			/** Define the constraints *con for later (repeated) solving. */
			virtual void set_constraints(const lin_constraint_system<
					scalar_type>& con) const {throw std::runtime_error("missing implementation");};

			/** Solve the linear program given by the objective function f under constraints
			 * defined previously with set_constraints(con). Formaly, compute max f(x) s.t. con(x).
			 * If f==0, check only satisfiability (support_vec,is_bounded,is_maximum undefined).
			 */
			virtual void maximize(const lin_expression<scalar_type>& f, typename lp_solver<scalar_type>::lp_result& res) const {throw std::runtime_error("missing implementation"); };


};

			}

#endif /*LP_SOLVER_MONNIAUX_H_*/
