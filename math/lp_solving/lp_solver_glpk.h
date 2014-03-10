#ifndef LP_SOLVER_GLPK_H_
#define LP_SOLVER_GLPK_H_

#include "utility/stl_helper_functions.h"
#include "math/lp_solving/lp_solver.h"
#include <glpk.h>
#include "glpk_cleaner.h"

namespace math {

template<typename scalar_type> class lp_solver_glpk: public lp_solver<scalar_type> {
public:
	typedef boost::shared_ptr<lp_solver_glpk<scalar_type> > ptr;
	typedef boost::shared_ptr<const lp_solver_glpk<scalar_type> > const_ptr;

	typename lp_solver<scalar_type>::lp_result lp_result;

	lp_solver_glpk() :
		my_glp_prob(0), my_constraints_empty(false) {
//		LOGGER_OS(DEBUG7, __FUNCTION__) << "created glpk solver " << this;
	}
	;

	/** Virtual Destructor. */
	virtual ~lp_solver_glpk() {
//		LOGGER_OS(DEBUG7, __FUNCTION__) << "destructing glpk solver " << this;
		if (my_glp_prob) {
			glpk_cleaner::delete_glp_prob(my_glp_prob);
			my_glp_prob = 0;
		}
	}
	;

	/** Solve the linear program given by the objective function f under the constraints
	 * con. Formaly, compute max f(x) s.t. con(x).
	 * @todo check strict inequalities
	 */
	virtual void maximize(const lin_expression<scalar_type>& f, const lin_constraint_system<
			scalar_type>& cons, typename lp_solver<scalar_type>::lp_result& res) {
		set_constraints(cons);
		maximize(f, res);
	}

	/** Define the constraints *con for later (repeated) solving. */
	virtual void set_constraints(const lin_constraint_system<scalar_type>& passed_cons) {
		if (passed_cons.size() > 0) {
			my_constraints_empty = false;

			/* If this is the first time, set parameters */
			if (!my_glp_prob) {
				// define the parameters
				glp_init_smcp(&my_glp_parm);
				my_glp_parm.msg_lev = GLP_MSG_OFF;
	// The following lead to GLP_UNDEF undefined solution on a perfectly fine instance.
	//				my_glp_parm.presolve = GLP_ON;
				my_glp_parm.tol_bnd = 1e-9;
			}
			create_glpk_problem();

			// if necessary, unify the domains
			const lin_constraint_system<scalar_type>* pcons=&passed_cons;
			lin_constraint_system<scalar_type> unified_cons;
			if (!passed_cons.has_unified_domains()) {
				unified_cons=passed_cons;
				unified_cons.unify_domains();
				pcons=&unified_cons;
			}
			const lin_constraint_system<scalar_type>& cons(*pcons);

//std::cout << cons << std::endl;
			positional_vdomain dom=cons.begin()->get_l().domain();

			unsigned int n=dom.size();
			unsigned int m=cons.size();

			int coeff_count = m*n;

			/* Get an iimap for the whole thing. */
			variable_id_set all_vars = cons.get_variable_ids();
			assert(set_contains(all_vars,dom.get_variable_ids()));
			assert(set_contains(dom.get_variable_ids(),all_vars));
			my_is_satisfiable = true;
			my_without_vars = false;
			if (all_vars.empty()) {
				my_without_vars = true;
				my_is_satisfiable = math::maybe(cons.is_satisfiable());
				return;
			}
			my_iimap = dom.get_index_to_variable_id_map();
//std::cout << my_iimap << std::endl;

			int ia[1 + coeff_count], ja[1 + coeff_count];
			double ar[1 + coeff_count];

			if (suppress_glpk_terminal_output)
				glp_term_out(GLP_OFF);
			else
				glp_term_out(GLP_ON);

			//glp_set_prob_name(my_glp_prob, "sspaceex");
			glp_set_obj_dir(my_glp_prob, GLP_MAX);
			glp_add_rows(my_glp_prob, cons.size());
			glp_add_cols(my_glp_prob, my_iimap->dimensions());
			/* Define all variables as free. */
			for (unsigned int i = 0; i < all_vars.size(); ++i) {
				glp_set_col_bnds(my_glp_prob, i + 1, GLP_FR, 0.0, 0.0);
			}

			/* Define the coefficients and bounds for each constraint.
			 * A constraint a1 x1 + ... + an xn + b SIGN 0 corresponds according
			 * to SIGN to the following bounds:
			 * - LT/LE : GLP_UP, 0.0,  -b
			 * - GT/GE : GLP_LO,  -b, 0.0
			 * - EQ    : GLP_FX,  -b, 0.0
			 * */
			unsigned int k=0;

			typename lin_constraint_system<scalar_type>::const_iterator it=cons.begin();
			for (unsigned int i=1;i<=m;++i) {
				const lin_expression<scalar_type>& l=it->get_l();

				for (unsigned int j=0;j<n;) {
					++k;
					ia[k] = i;
					ar[k] = convert_element<double>(l[j]);
					++j;
					ja[k] = j;
				}

				if (is_LT_or_LE(it->get_sign())) {
					glp_set_row_bnds(my_glp_prob, i, GLP_UP, 0.0,
							-double(it->get_inh_coeff()));
				} else if (is_GT_or_GE(it->get_sign())) {
					glp_set_row_bnds(my_glp_prob, i, GLP_LO,
							-double(it->get_inh_coeff()), 0.0);
				} else {
					// EQ
					glp_set_row_bnds(my_glp_prob, i, GLP_FX,
							-double(it->get_inh_coeff()), 0.0);
				}

				++it;
			}
			assert(k==coeff_count);

			/* Set the constraints */
			glp_load_matrix(my_glp_prob, coeff_count, ia, ja, ar);
			//glp_write_lp(my_glp_prob, NULL, "/dev/stdout");
		} else {
			my_constraints_empty = true;
			my_iimap = index_to_variable_id_map_ptr();
			my_is_satisfiable = true;
			my_without_vars = true;
		}
	}
	;

	/** Solve the linear program given by the objective function f under constraints
	 * defined previously with set_constraints(con). Formaly, compute max f(x) s.t. con(x).
	 * If f==0, check only satisfiability (support_vec,is_bounded,is_maximum undefined).
	 */
	virtual void maximize(const lin_expression<scalar_type>& f_orig,
			typename lp_solver<scalar_type>::lp_result& res) {
		LOGGERSWOC(DEBUG4,"lp_solver_glpk::maximize","Solving LP with glpk");

		if (my_constraints_empty) {
			/** No constraints. */
			res.is_bounded = f_orig.is_zero();
			res.is_maximum = false;
			res.is_unsat = false;
		} else if (my_without_vars) {
			if (my_is_satisfiable) {
				res.is_bounded = f_orig.is_zero();
				res.is_maximum = false;
				res.is_unsat = false;
			} else {
				res.is_bounded = true;
				res.is_maximum = false;
				res.is_unsat = true;
			}
		} else {
			if (my_glp_prob) {
				bool f_has_unconstrained_vars=false;
//std::cout << "cost:" << f_orig << std::endl;
//std::cout << "with iimap:" << f_orig.get_index_to_variable_id_map() << std::endl;
//std::cout << "my_iimap:" << my_iimap << std::endl;
				lin_expression<scalar_type> f=f_orig;
				if(!my_iimap->contains_variables(f.get_index_to_variable_id_map())) {
					// if there are variables that are unconstrained (not in my_iimap)
					// and they have non-zero coeffs, the problem is unbounded
					// if they have zero coeffs, remove them

					for (unsigned int i = 0; i < f_orig.size() && !f_has_unconstrained_vars; i++) {
						if (f_orig[i] != scalar_type(0) && !my_iimap->has_id(f_orig.get_id(i))) {
							f_has_unconstrained_vars=true;
						}
					}
					f.remap(my_iimap);
				}

				/* Reset objective function to zero */
				for (unsigned int i = 0; i < my_iimap->dimensions(); ++i) {
					glp_set_obj_coef(my_glp_prob, i + 1, 0.0);
				}

				/* Set the objective function to passed linear expression */
				if (my_iimap == f.get_index_to_variable_id_map()) {
					for (unsigned int i = 0; i < f.size(); ++i) {
//						std::cout << f.get_variable(i) << "=" << double(f[i]) << ",";
						glp_set_obj_coef(my_glp_prob, i + 1, double(f[i]));
					}
				} else {
					for (unsigned int i = 0; i < f.size(); ++i) {
						glp_set_obj_coef(my_glp_prob, my_iimap->get_index(f.get_id(i)) + 1,
								double(f[i]));
					}
				}

				glp_set_obj_coef(my_glp_prob, 0, double(f.get_inh_coeff()));

				//glp_write_lp(my_glp_prob, NULL, "/dev/stdout");

				/* carry out the optimization */
				glp_simplex(my_glp_prob, &my_glp_parm);

				/* get the optimal objective function value */
				res.max_value = scalar_type(glp_get_obj_val(my_glp_prob));

				/* get the optimal point */
				res.support_vec = vdom_vector<scalar_type> (my_iimap);
				for (unsigned int i = 0; i < my_iimap->dimensions(); ++i) {
					res.support_vec[i] = scalar_type(glp_get_col_prim(my_glp_prob, i + 1));
				};

//glp_write_lp(my_glp_prob, NULL, "/dev/stdout");

				/* get info on the problem */
				int glp_stat = glp_get_status(my_glp_prob);
				/*
				 GLP_OPT    solution is optimal;
				 GLP_FEAS   solution is feasible;
				 GLP_INFEAS solution is infeasible;
				 GLP_NOFEAS problem has no feasible solution;
				 GLP_UNBND  problem has unbounded solution;
				 GLP_UNDEF  solution is undefined.
				 */
				res.is_unsat = false;
				res.is_bounded = true;
				res.is_maximum = true;
				if (glp_stat == GLP_OPT) { // everything is fine
				} else if (glp_stat == GLP_NOFEAS) {
					res.is_unsat = true;
				} else if (glp_stat == GLP_UNBND) {
					res.is_bounded = false;
					res.is_maximum = false;
				} else {
					std::string resp="unknown";
					if (glp_stat == GLP_FEAS) { resp="GLP_FEAS";
					} else if (glp_stat == GLP_INFEAS) { resp="GLP_INFEAS";
					} else if (glp_stat == GLP_UNDEF) { resp="GLP_UNDEF";
					}
					glp_write_lp(my_glp_prob,NULL,"/dev/stdout");
					throw std::runtime_error("glpk " + std::string(glp_version()) + "returned status "
							+ int2string(glp_stat) + "("+resp+").");
				}
//				std::cout << std::endl << "GLPK LP: ";
//				glp_write_lp(my_glp_prob,NULL,"/dev/stdout");

				// If there are unconstrained variables in the objective function,
				// set the problem as unbounded, unless it is unsat
				if (!res.is_unsat && f_has_unconstrained_vars) {
					res.is_bounded = false;
					res.is_maximum = false;
				}

				/* after solving, set presolving off so the next time it won't
				 * be done unless new constraints are defined. */
				my_glp_parm.presolve = GLP_OFF;
			} else {
				throw std::runtime_error("lp constraints are undefined");
			}
		}
	}
	;

	/** Get a string name for the solver. */
	virtual std::string get_solver_name() const {
		return "lp_solver_glpk(v" + std::string(glp_version()) + ")";
	}
	;

	static bool suppress_glpk_terminal_output;

private:
	/** Create a glpk problem */
	void create_glpk_problem() {
		/* If there is already a problem defined, delete it first. */
		if (my_glp_prob) {
			glpk_cleaner::delete_glp_prob(my_glp_prob);
			my_glp_prob = 0;
		}
		my_glp_prob = glpk_cleaner::create_glp_prob();
//		LOGGER_OS(DEBUG7, __FUNCTION__) << "created glp problem " << my_glp_prob << " for solver " << this;
	}

	/** Create a glpk problem */
	void destroy_glpk_problem() {
		/* If there is already a problem defined, delete it first. */
		if (my_glp_prob) {
//			LOGGER_OS(DEBUG7, __FUNCTION__) << "destroying glp problem " << my_glp_prob << " of solver " << this;
			glpk_cleaner::delete_glp_prob(my_glp_prob);
			my_glp_prob = 0;
		}
	}

	glp_prob* my_glp_prob;
	glp_smcp my_glp_parm;
	bool my_constraints_empty;
	index_to_variable_id_map_ptr my_iimap;
	bool my_is_satisfiable;
	bool my_without_vars;
};

template<typename T> bool lp_solver_glpk<T>::suppress_glpk_terminal_output = true;

}

#endif /*lp_solver_glpk_H_*/
