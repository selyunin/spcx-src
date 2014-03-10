#ifndef FM_ELIMINATION_H_
#define FM_ELIMINATION_H_

#include <list>
#include "math/vdom/lin_constraint.h"

namespace math {
namespace fourier_motzkin {

/* Combine the constraints con1 and con2 such that the resulting coefficient of variable vid is zero.
 * The sign of the resulting constraint is the same direction as that of con1, i.e.,
 * con1 is always multiplied with a positive scalar.
 */
template<typename scalar_type> lin_constraint<scalar_type> eliminate_coeff(
		const lin_constraint<scalar_type>& con1,
		const lin_constraint<scalar_type>& con2, variable_id vid) {
	const scalar_type& c1=con1.get_coeff_with_id(vid);
	const scalar_type& c2=con2.get_coeff_with_id(vid);
	if (c1==scalar_type(0))
		return con1;
	else if (c2==scalar_type(0))
		return con2;
	else {
		lin_constraint<scalar_type> con;
		if (c2>=scalar_type(0))
			con=c2*con1+(-c1)*con2;
		else
			con=(-c2)*con1+c1*con2;
		//std::cout << con1 << " and " << con2 << " give " << con;
		con.set_coeff_with_id(vid, scalar_type(0)); // Not formally necessary, but to avoid numerical errors
		//std::cout << " zeroed " << con << std::endl;
		return con;
	}
}
;

/** Eliminate trivially redundant constraints from the list. */
template<typename scalar_type, typename passed_list_type> void eliminate_trivial_constraints(
		passed_list_type& con_sys) {
	typedef passed_list_type constraint_list_type;

	bool is_empty=false;
	typename constraint_list_type::iterator it=con_sys.begin();
	while (it!=con_sys.end() && !is_empty) {
		if (math::maybe(it->is_always_satisfied())) // the constraint is redundant
		{
			it=con_sys.erase(it);
		} else {
			is_empty=is_empty || !math::maybe(it->is_satisfiable());
			++it;
		}
	}
	if (is_empty) // replace by a single canonical constraint
	{
		con_sys=constraint_list_type();
		con_sys.push_back(lin_constraint<scalar_type>::unsatisfiable_constraint());
	}
}
;

/** Use Fourier-Motzkin elimination to eliminate the variable with id vid from the set of
 * constraints con_sys.
 * passed_list_type is the type of the set of linear constraints. It needs to provide
 * const_iterator, begin(), end() that dereference to lin_constraint<scalar_type> with
 * -> and *.
 *
 * @attention For efficiency reasons, the variables in vid are not removed from the
 * index_to_variable_id_maps of the constraints.
 * */
template<typename scalar_type, typename passed_list_type> passed_list_type eliminate(
		const passed_list_type& con_sys, variable_id vid) {
	/* optionally : if scalar_type is approximative, make sure that
	 * very small coefficients are replaced by zeros. */
	//	typedef typename std::list<lin_constraint<scalar_type> >
	typedef passed_list_type constraint_list_type;
	constraint_list_type new_cons;

	/* Collect constraints whose coefficient of vid is positive, negative, or zero. */
	constraint_list_type pos_cons, neg_cons, zero_cons, eq_cons;
	for (typename constraint_list_type::const_iterator it=con_sys.begin(); it
			!=con_sys.end(); ++it) {
		const scalar_type& a=it->get_canonic_coeff_with_id(vid);
		if (a==scalar_type(0))
			zero_cons.push_back(*it);
		else if (it->is_inequality()) {
			if (a>scalar_type(0))
				pos_cons.push_back(*it);
			else
				// if (a<scalar_type(0))
				neg_cons.push_back(*it);
		} else
			eq_cons.push_back(*it);
	}
	//std::cout << pos_cons.size() << neg_cons.size() << zero_cons.size() << eq_cons.size();

	/* Combine positive and negative inequality constraints. */
	if (pos_cons.size()>0 && neg_cons.size()>0) {
		for (typename constraint_list_type::const_iterator jt=pos_cons.begin(); jt
				!=pos_cons.end(); ++jt) {
			for (typename constraint_list_type::const_iterator kt=
					neg_cons.begin(); kt !=neg_cons.end(); ++kt) {
				new_cons.push_back(eliminate_coeff(*jt, *kt, vid));
			}
		}
	}
	/* Combine equalities with positive and negative inequality constraints. */
	for (typename constraint_list_type::const_iterator it=eq_cons.begin(); it
			!=eq_cons.end(); ++it) {
		for (typename constraint_list_type::const_iterator jt= pos_cons.begin(); jt
				!=pos_cons.end(); ++jt) {
			new_cons.push_back(eliminate_coeff(*jt, *it, vid));
		}
		for (typename constraint_list_type::const_iterator jt= neg_cons.begin(); jt
				!=neg_cons.end(); ++jt) {
			new_cons.push_back(eliminate_coeff(*jt, *it, vid));
		}
	}
	/* Combine equalities with the other equalities. . */
	for (typename constraint_list_type::const_iterator it=eq_cons.begin(); it
			!=eq_cons.end(); ++it) {
		typename constraint_list_type::const_iterator jt=it;
		++jt;
		while (jt!=eq_cons.end()) {
			new_cons.push_back(eliminate_coeff(*jt, *it, vid));
			++jt;
		}
	}

	//std::cout << "new " << new_cons << std::endl;
	/* Eliminate any trivial constraints in case they were introduced. */
	eliminate_trivial_constraints<scalar_type>(new_cons);
	//std::cout << "red elim " << new_cons << std::endl;

	/* Add constraints for which the coefficient was zero. */
	for (typename constraint_list_type::const_iterator it=zero_cons.begin(); it
			!=zero_cons.end(); ++it) {
		new_cons.push_back(*it);
	}
	return new_cons;
}
;

template<typename parameter_list_type> variable_id_set get_variable_ids(
		const parameter_list_type& con_sys) {
	/* to find a variable_id that's not yet in con_sys, get all variables */
	variable_id_set vis;
	for (typename parameter_list_type::const_iterator it=con_sys.begin(); it
			!=con_sys.end(); ++it) {
		variable_id_set con_vars=(*it).get_variable_ids();
		vis.insert(con_vars.begin(), con_vars.end());
	}
	return vis;
}
;

/** Use elimination to find the bounds on linear cost function.
 * Note that this function is extremely slow and only for experimental purposes.
 *
 * parameter_list_type is the type of the set of linear constraints. It needs to provide
 * const_iterator, begin(), end() that dereference to lin_constraint<scalar_type> with
 * -> and *.
 * 
 * If the parameter list is empty, then set lower_bounded = upper_bound = is_empty = false.
 */
template<typename scalar_type, typename parameter_list_type> 
	void 
	optimize(const parameter_list_type& con_sys, 
	         lin_expression<scalar_type> f,
	         scalar_type& min_val, scalar_type& max_val, 
	         bool& lower_bounded,	bool& upper_bounded, 
	         bool& lower_strict,  bool& upper_strict,
	         bool& is_empty) {

	typedef typename std::list<lin_constraint<scalar_type> >
			constraint_list_type;


	/*
	  check empty first
	 */
	if (con_sys.begin() == con_sys.end()) {
		lower_bounded = false;
		upper_bounded = false;
		is_empty = false;
		return;
	}
	/* to find a variable_id that's not yet in con_sys, get all variables */
	variable_id_set vis = get_variable_ids(con_sys);

	variable_id cost_id = variable::get_temp_variable_id();
	constraint_list_type cost_cons;
	for (typename parameter_list_type::const_iterator it=con_sys.begin(); it
			!=con_sys.end(); ++it) {
		cost_cons.push_back(*it);
	}
	/* Add the constraint cost_id==f, i.e., f+(-1)*cost_id==0 */
	lin_constraint<scalar_type> cost_con(f, EQ);
	cost_con.set_coeff_with_id(cost_id, scalar_type(-1));
	cost_cons.push_back(cost_con);

	//std::cout << "before elim : " << cost_cons << std::endl;
	/* Eliminate all variables. */
	for (variable_id_set::const_iterator vit = vis.begin(); vit != vis.end(); ++vit) {
		cost_cons=eliminate<scalar_type, constraint_list_type>(cost_cons, *vit);
		//std::cout << "after elim " << *vit << ": " << cost_cons << std::endl;
	}

	/* parse the constraints to find upper and lower bounds */
	lower_bounded=false;
	upper_bounded=false;
	lower_strict=false;
	upper_strict=false;
	is_empty=false;
	for (typename constraint_list_type::const_iterator it = cost_cons.begin(); it
			!= cost_cons.end() && !is_empty; ++it) {
		//std::cout << std::endl << *it << std::endl;
		scalar_type c=it->get_canonic_coeff_with_id(cost_id); // c*x+b<=0
		scalar_type b=it->get_canonic_inh_coeff();
		if (c!=scalar_type(0)) {
			if (it->is_equality() || c>scalar_type(0)) // upper bound
			{
				if (!upper_bounded || -b/c<=max_val) {
					max_val=-b/c;
					upper_bounded=true;
					upper_strict=upper_strict || it->is_strict_inequality();
					//std::cout << max_val << " upper" << upper_strict << std::endl;
				}
			} else if (it->is_equality() || c<scalar_type(0)) // lower bound
			{
				if (!lower_bounded || -b/c>=min_val) {
					min_val=-b/c;
					lower_bounded=true;
					lower_strict=lower_strict || it->is_strict_inequality();
					//std::cout << min_val << " lower" << lower_strict << std::endl;
				}
			}
		} else // c==0
		{
			is_empty=is_empty || !math::maybe(it->is_satisfiable());
		}
		//std::cout << is_empty << std::endl;
	}
	is_empty=is_empty || (lower_bounded && upper_bounded && (max_val<min_val
			|| (max_val==min_val && (lower_strict || upper_strict))));

	/* remove the temp variable */
	variable::remove_variable_id(cost_id);
}
;

}
}

#endif /*FM_ELIMINATION_H_*/
