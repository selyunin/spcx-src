/*
 * Algebraic Variable Elimination
 *
 *  Created on: Dec 6, 2010
 *      Author: Ray
 */

#ifndef ALG_VAR_ELIM_H_
#define ALG_VAR_ELIM_H_

#include "math/vdom/lin_constraint_system.h"
#include "math/vdom/lin_constraint_system_utility.h"
#include "math/matrix.h"

namespace math {

/**
 * Class having the routine to eliminate required set of variables (u_id_set and
 * v_id_set) from a given set of linear constraints.
 */

template<typename scalar_type>
class alg_var_elim{

public:
	alg_var_elim(positional_vdomain xdom, positional_vdomain udom, positional_vdomain vdom,
			matrix<scalar_type> A_mat, matrix<scalar_type> B_mat, vector<scalar_type> alpha_vec,
			const lin_constraint_system<scalar_type>& cons){
		assert(A_mat.size1() == A_mat.size2());
		assert(A_mat.size1() == B_mat.size1());
		assert(A_mat.size1() == alpha_vec.size());
		assert(A_mat.size2() == xdom.size());
		assert(B_mat.size2() == udom.size());

		n = xdom.size();
		x_dom = xdom;
		u_dom = udom;
		v_dom = vdom;
//		std::cout << "xdom:" << x_dom << " udom:" << u_dom << " vdom:" << v_dom << std::endl;
		A = A_mat;
		B = B_mat;
		alpha = alpha_vec;
		all_cons = cons;
		all_cons.collapse_inequalities();
		nb_cons = all_cons.size();
		nb_ineq = 0;
		assert(C.size1() == D.size1() && D.size1() == E.size1() && E.size1() == beta.size());
	};
	/**
	 * Attempts to eliminate the algebraic variables.
	 * Returns A_elim the matrix after the elimination algorithm
	 *         is applied to the matrix of the coefficients of the system variables
	 *         B_elim the matrix after the elimination algorithm
	 *         is applied to the matrix of the coefficients of the algebraic variables.
	 *         B_elim is a zero matrix in the best possible case.
	 *         reduced_cons are the constraints on the variables (after reduction if possible)
	 */
	void eliminate(matrix<scalar_type>& A_elim, matrix<scalar_type>& B_elim, positional_vdomain& B_dom, vector<scalar_type>& alpha_elim, lin_constraint_system<scalar_type>& reduced_cons) ;


	/*
	 * Returns the matrices from the coefficients of the invariant constraint.
	 * C is the matrix of the coefficients of the system variables
	 * D is the matrix of the coefficients of the algebraic variables
	 * E is the matrix of the coefficients of the other variables
	 */
	void get_cons_matrices(matrix<scalar_type>&C, matrix<scalar_type>&D,matrix<scalar_type>&E) const ;
	void print() const;

private:
	unsigned int n,nb_cons,nb_ineq;

	matrix<scalar_type> A;
	matrix<scalar_type> B;
	matrix<scalar_type> C;
	matrix<scalar_type> D;
	matrix<scalar_type> E;
	vector<scalar_type> alpha;
	vector<scalar_type> beta;
	typedef typename lin_constraint<scalar_type>::sign sign;
	std::vector<sign> ineq_signs;

	matrix<scalar_type> A_elim; // system vars coeff matrix after elimination
	matrix<scalar_type> B_elim; // algebraic vars coeff matrix after elimination
	vector<scalar_type> alpha_elim; // inh coeffiecients after elimination process
	positional_vdomain x_dom; // system variables
	positional_vdomain u_dom; // algebraic variables
	positional_vdomain v_dom; // other variables
	positional_vdomain slack_dom; // slack variables
	positional_vdomain vs_dom; // other+slack variables

	lin_constraint_system<scalar_type> disregarded_cons;
	lin_constraint_system<scalar_type> all_cons;

	/** Retrieve matrix equalities */
	void get_matrix_form_from_equalities(lin_constraint_system<scalar_type> cons);

	/** Retrieve matrices by introducing slack variables
	 */
	void get_matrix_form_with_slack(lin_constraint_system<scalar_type> cons);
};

}

#include "alg_var_elim.hpp"


#endif
