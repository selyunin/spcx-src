/*
 * lin_constraint_system_utility.h
 *
 *  Created on: May 8, 2010
 *      Author: frehse
 */

#ifndef LIN_CONSTRAINT_SYSTEM_UTILITY_H_
#define LIN_CONSTRAINT_SYSTEM_UTILITY_H_

#include "lin_constraint_system.h"
#include "math/matrix.h"

namespace math {

/** Return the system in canonic matrix form A*x<=b with domain dom.
 *
 * Equality constraints are converted to inequality constraints.
 * Strict inequality constraints are treated as nonstrict.
 * */
template<typename scalar_type>
void canonic_matrix_form(matrix<scalar_type>& A, vector<scalar_type>& b,
		positional_vdomain& dom, lin_constraint_system<scalar_type> cons) {

	if (cons.size()==0) {
		A=matrix<scalar_type>();
		b=vector<scalar_type>();
		dom=positional_vdomain();
	} else {
		cons.unify_domains();
		cons.expand_equalities();
		dom=cons.begin()->get_l().domain();

		unsigned int n=dom.size();
		unsigned int m=cons.size();

		A=matrix<scalar_type>(m,n);
		b=vector<scalar_type>(m);

		typename lin_constraint_system<scalar_type>::const_iterator it=cons.begin();
		for (unsigned int i=0;i<m;++i) {
			lin_expression<scalar_type> l=it->get_canonic_l();
			for (unsigned int j=0;j<n;++j) {
				A(i,j)=l[j];
			}
			b[i]=-it->get_canonic_inh_coeff();
			++it;
		}
	}
}
;

/** Construct a system of linear constraints from matrix form
 *
 * Constructs the system A*x<=b given A,b and a domain.
 * Passed a nonempty vector of signs, the default sign <= is replaced
 * for with the corresponding sign.
 */
template<typename scalar_type>
lin_constraint_system<scalar_type> construct_linear_constraint_system(
		const matrix<scalar_type>& A, const vector<scalar_type>& b,
		const positional_vdomain& dom,
		const std::vector<typename lin_constraint<scalar_type>::sign>& signs =
				std::vector<typename lin_constraint<scalar_type>::sign>()) {
	assert(A.size1() == b.size());
	assert(A.size2() == dom.size());
	assert(signs.empty() || A.size1() == signs.size());

	bool use_signs = !signs.empty();

	lin_constraint_system<scalar_type> cons;
	typename lin_constraint<scalar_type>::sign s;
	for (size_t i = 0; i < A.size1(); ++i) {
		vdom_vector<scalar_type> vec(dom, A.vector_from_row(i));
		if (use_signs) {
			s = signs[i];
		} else {
			s = LE;
		}
		lin_constraint<scalar_type> con(vec, -b[i], s);
		if (!maybe(con.is_always_satisfied())) {
			cons.insert(con);
		}
	}
	return cons;
}

/** Computes the variables with nonzero coefficients in the constraints */
template<typename scalar_type>
variable_id_set compute_used_variables(const lin_constraint_system<scalar_type>& cons) {
	variable_id_set vis;
	for (typename lin_constraint_system<scalar_type>::const_iterator it = cons.begin(); it!=cons.end(); ++it) {
		variable_id_set local_vis = it->get_l().get_used_variable_ids();
		set_union_assign(vis,local_vis);
	}
	return vis;
}

/** Removes unused variables from the domains of the constraints */
template<typename scalar_type>
void remove_unused_variables(lin_constraint_system<scalar_type>& cons) {
	variable_id_set used = compute_used_variables(cons);
	variable_id_set vars = cons.get_variable_ids();
	set_difference_assign(vars,used);
	cons.remove_variables(vars);
}

}

#endif /* LIN_CONSTRAINT_SYSTEM_UTILITY_H_ */
