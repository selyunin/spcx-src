#ifndef LP_SOLVER_FM_H_
#define LP_SOLVER_FM_H_

#include "utility/stl_helper_functions.h"
#include "math/lp_solving/lp_solver.h"
#include "math/lp_solving/fourier_motzkin/fm_elimination.h"

namespace math {

template<typename scalar_type> class lp_solver_fm: public lp_solver<scalar_type> {
public:
	typedef boost::shared_ptr<lp_solver_fm<scalar_type> > ptr;
	typedef boost::shared_ptr<const lp_solver_fm<scalar_type> > const_ptr;

	typename lp_solver<scalar_type>::lp_result lp_result;

	lp_solver_fm() :
		my_constraints(NULL) {
	}
	;

	/** Virtual Destructor. */
	virtual ~lp_solver_fm() {
	}
	;

	/** Solve the linear program given by the objective function f under the constraints
	 * con. Formaly, compute max f(x) s.t. con(x).
	 */
	virtual void maximize(const lin_expression<scalar_type>& f, const lin_constraint_system<
			scalar_type>& con_sys, typename lp_solver<scalar_type>::lp_result& res) {
		scalar_type min_val; // ignored
		bool lower_bounded; // ignored
		bool lower_strict; // ignored
		bool upper_strict; // negated
		fourier_motzkin::optimize(con_sys, f, min_val, res.max_value, lower_bounded,
				res.is_bounded, lower_strict, upper_strict, res.is_unsat);
		res.is_maximum = !upper_strict;
		res.support_vec.clear();
		//std::cout << "max " << f << " s.t. " << con_sys << "="; print<scalar_type>(std::cout,res);
		if (!f.is_zero()) {
			/* Find a support vector by fixing one variable at a time */
			variable_id_set vis = fourier_motzkin::get_variable_ids(con_sys);

			/* Get the set of constraints, take the closure */
typedef			typename std::list<lin_constraint<scalar_type> >
			constraint_list_type;
			constraint_list_type cons;
			for (typename lin_constraint_system<scalar_type>::const_iterator it=
					con_sys.begin(); it !=con_sys.end(); ++it) {
				lin_constraint<scalar_type> con=*it;
				con.closure_assign();
				cons.push_back(con);
			}
			//std::cout << "f " << f << " max " << res.max_value << std::endl;
			/* Add the optimal cost value as a constraint,
			 * i.e., f==max_val <=>  f-max_val==0. */
			lin_expression<scalar_type> l=f;
			l.set_inh_coeff(f.get_inh_coeff()-res.max_value);
			cons.push_back(lin_constraint<scalar_type>(l, EQ));

			variable_id_set::const_iterator vit=vis.begin();
			while (vit!=vis.end()) {
				//std::cout << "var " << *vit << " cons: " << cons << std::endl;
				lin_expression<scalar_type> f2;
				f2.set_coeff_with_id(*vit, scalar_type(1));
				scalar_type max_val;
				bool upper_bounded, is_empty;
				fourier_motzkin::optimize(cons, f2, min_val, max_val,
						lower_bounded, upper_bounded, lower_strict,
						upper_strict, is_empty);
				/* fix the value of variable *vit to max_val */
				res.support_vec.set_coeff_with_id(*vit, max_val);
				/* add the constraint *vit == max_val */
				lin_expression<scalar_type> l2;
				l2.set_coeff_with_id(*vit, scalar_type(1));
				l2.set_inh_coeff(-max_val);
				cons.push_back(lin_constraint<scalar_type>(l, EQ));
				/* go on to the next variable */
				++vit;
			}
		}
	}
	;

	/** Define the constraints *con for later (repeated) solving. */
	virtual void set_constraints(const lin_constraint_system<
			scalar_type>& con) {my_constraints=&con;};

	/** Solve the linear program given by the objective function f under constraints
	 * defined previously with set_constraints(con). Formaly, compute max f(x) s.t. con(x).
	 * If f==0, check only satisfiability (support_vec,is_bounded,is_maximum undefined).
	 */
	virtual void maximize(const lin_expression<scalar_type>& f, typename lp_solver<scalar_type>::lp_result& res) {
		if (my_constraints) {
			maximize(f,*my_constraints,res);
		}
		else
		{
			throw std::runtime_error("constraints not defined");
		}
	};

	/** Get a string name for the solver. */
	virtual std::string get_solver_name() const {
		return "lp_solver_fm";
	};

private:
	const lin_constraint_system<scalar_type>* my_constraints;
};

}

#endif /*lp_solver_fm_H_*/
