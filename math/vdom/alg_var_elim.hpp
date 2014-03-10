
namespace math {

template<typename scalar_type>
void alg_var_elim<scalar_type>::get_matrix_form_from_equalities(lin_constraint_system<scalar_type> cons){

	disregarded_cons = lin_constraint_system<scalar_type>();
	ineq_signs = std::vector<sign>();
	positional_vdomain dom;

	if (cons.size()==0) {
		C = matrix<scalar_type>();
		D = matrix<scalar_type>();
		E = matrix<scalar_type>();
		dom = positional_vdomain();
		nb_cons = 0; // # rows
		nb_ineq = 0;
	} else {
		dom = cons.unify_domains();
		nb_cons = cons.size(); // # rows

		// count the number of equalities (we don't consider inequalities for now)
		size_t nb_eq = 0;
		for (typename lin_constraint_system<scalar_type>::const_iterator it=cons.begin();it!=cons.end();++it) {
			if (it->is_equality()) {
				++nb_eq;
			}
		}
		nb_ineq = nb_cons-nb_eq; // # inequalities

		beta = vector<scalar_type>(nb_cons);

		C = matrix<scalar_type>(nb_cons, x_dom.size());
		D = matrix<scalar_type>(nb_cons, u_dom.size());
		E = matrix<scalar_type>(nb_cons, v_dom.size());

		// Keep track of the constraints we don't take into account
		typename lin_constraint_system<scalar_type>::const_iterator it=cons.begin();
		lin_expression<scalar_type> lx,lu,lv;
		size_t i = 0;
		for (; it != cons.end(); ++it) {
			if (it->is_equality()) {
				lin_expression<scalar_type> l = it->get_canonic_l();
				lx = l;
				lu = l;
				lv = l;
				lx.remap(x_dom);
				C.assign_row(i, lx.get_vector());
				lu.remap(u_dom);
				D.assign_row(i, lu.get_vector());
				lv.remap(v_dom);
				E.assign_row(i, lv.get_vector());
				beta[i] = it->get_canonic_inh_coeff();
				++i;
			} else {
				disregarded_cons.insert(*it);
			}
		}

	}
}

template<typename scalar_type>
void alg_var_elim<scalar_type>::get_matrix_form_with_slack(lin_constraint_system<scalar_type> cons){

	// there are no disregarded constraints
	disregarded_cons = lin_constraint_system<scalar_type>();
	ineq_signs = std::vector<sign>();
	positional_vdomain dom;

	dom = cons.unify_domains();
	nb_cons = cons.size();
	slack_dom = positional_vdomain();
	if (cons.size()==0) {
		C = matrix<scalar_type>();
		D = matrix<scalar_type>();
		E = matrix<scalar_type>();
		nb_ineq = 0;
		beta = vector<scalar_type>();
	} else {

		// count the number of equalities (we don't consider inequalities for now)
		size_t nb_eq = 0;
		for (typename lin_constraint_system<scalar_type>::const_iterator it=cons.begin();it!=cons.end();++it) {
			if (it->is_equality()) {
				++nb_eq;
			}
		}
		nb_ineq = nb_cons-nb_eq; // # inequalities

		// @todo Check if there are already slack variables

		// Create slack variables
		for (size_t i = 0; i< nb_ineq; ++i) {
			slack_dom.add_variable(variable("SLACK"+to_string(i)));
		}
		LOGGER_OS(DEBUG7,__FUNCTION__) << "slack domain: " << slack_dom << std::endl;
		// Add slack variables to the domain of v
		vs_dom = compose(v_dom,slack_dom);

		C = matrix<scalar_type>(nb_cons, x_dom.size());
		D = matrix<scalar_type>(nb_cons, u_dom.size());
		E = matrix<scalar_type>(nb_cons, vs_dom.size());
		beta = vector<scalar_type>(nb_cons);
		ineq_signs = std::vector<sign>(nb_ineq);

		// Keep track of the constraints we don't take into account
		typename lin_constraint_system<scalar_type>::const_iterator it=cons.begin();
		lin_expression<scalar_type> lx,lu,lv;
		size_t i = 0;
		size_t ineq_count = 0;
		for (; it != cons.end(); ++it) {
			lin_expression<scalar_type> l = it->get_canonic_l();
			lx = l;
			lu = l;
			lv = l;
			lx.remap(x_dom);
			C.assign_row(i, lx.get_vector());
			lu.remap(u_dom);
			D.assign_row(i, lu.get_vector());
			lv.remap(vs_dom);
			if (it->is_inequality()) {
				// add slack variable
				const variable& slack_var = slack_dom.get_variable(ineq_count);
				lv.set_existing_coeff(slack_var,scalar_type(1));
				ineq_signs[ineq_count] = it->get_canonic_sign();
				++ineq_count;
			}
			E.assign_row(i, lv.get_vector());
			beta[i] = it->get_canonic_inh_coeff();

			++i;
		}
	}
	assert(slack_dom.size()==nb_ineq);
}
template<typename scalar_type>
void alg_var_elim<scalar_type>::eliminate(matrix<scalar_type>& A_elim, matrix<scalar_type>& B_elim, positional_vdomain& B_dom, vector<scalar_type>& alpha_elim,lin_constraint_system<scalar_type>& reduced_cons){
	using namespace numeric;

	// start with an empty set of constraints
	reduced_cons = lin_constraint_system<scalar_type>();

	// first try the simple method
	get_matrix_form_from_equalities(all_cons);

	bool success = false;
	bool D_is_singular;
	matrix<scalar_type> D_inv;
	// check if D is square and call inverse if only square
	if(D.size1() == D.size2())
		D_inv = D.inverse(D_is_singular);

	//check if E is 0 matrix and D is invertible
	if(E.is_zero() && !D_is_singular){
		LOGGER_OS(DEBUG6,__FUNCTION__) << "D invertible" << std::endl;
		A_elim = A - B*D_inv*C;
		// return 0 matrix.
		B_elim = matrix<scalar_type>(B.size1(),B.size2());
		alpha_elim = alpha - B * D_inv * beta;
		success = true;
	} else if(D.size1() == E.size1() && D.size2() + E.size2() == D.size1()){
		LOGGER_OS(DEBUG6,__FUNCTION__) << "DE invertible?" << std::endl;
		// The condition checks if the augmented matrix [D|E] is square
		matrix<scalar_type> DE = augment_cols(D,E);
		assert(DE.size1() == DE.size2());

		matrix<scalar_type> DE_inv;
		bool singular;
		DE_inv = DE.inverse(singular);
		if(!singular){
			LOGGER_OS(DEBUG6,__FUNCTION__) << "DE invertible!" << std::endl;
			matrix<scalar_type> zero(A.size1(), DE.size1() - B.size2());
			matrix<scalar_type> B0 = augment_cols(B,zero);

			matrix<scalar_type> phi = B0 * DE_inv;
			A_elim = A - (phi * C);
			B_elim = matrix<scalar_type>(B.size1(),B.size2());// zero matrix

			alpha_elim = alpha - phi * beta;

			// Missing: -(0 F)(D E)⁻1(Cx+d) >= 0
			// @todo Here we could return the reduced constraints for use as invariant
			success = true;
		}
	}
	if (!success) {
		LOGGER_OS(DEBUG6,__FUNCTION__) << "elimination with row echelon form" << std::endl;
		get_matrix_form_with_slack(all_cons);

		// Build the matrix to be reduced
		matrix<scalar_type> M(nb_cons,x_dom.size()+u_dom.size()+vs_dom.size()+1,scalar_type(0));

//		index_to_variable_id_map_ptr xmap = x_dom.get_index_to_variable_id_map();
//		positional_vdomain M_dom(xmap->get_map_with_primedness_increased());
		positional_vdomain w_dom = compose(u_dom,vs_dom);
		positional_vdomain M_dom = compose(w_dom,x_dom);
		assert(w_dom.size()==u_dom.size()+vs_dom.size());
		assert(M_dom.size()==w_dom.size()+x_dom.size());

		M_dom.add_variable(variable("CONST"));

		// Assume the equations have the form
		//		x' = Ax + Bu + b
		//		0  = Cx + Du + Ev + d
		//      0 <= Fv
		// here: alpha = b, beta = d
		matrix<scalar_type> F=diagonal_matrix(n,scalar_type(1));
		// Eliminate from the following form:
		//   u v x 1
		// [ D E C d]
		size_t i,j,u,v;
		vector<scalar_type> last_column(nb_cons);
		// first row
		i=0;
		j=0;
		u=D.size2();
		M.submatrix_assign(D,j,j+D.size1(),i,i+u);
		i+=u;
		u=E.size2();
		M.submatrix_assign(E,j,j+E.size1(),i,i+u);
		i+=u;
		u=C.size2();
		M.submatrix_assign(C,j,j+C.size1(),i,i+u);
		i+=u;
		last_column.subvector_assign(beta,j,j+beta.size());
		M.assign_column(i,last_column);

//std::cout << M << std::endl;
		reduced_row_echelon_form(M);
//std::cout << "reduced: " << std::endl << M << std::endl;

		// Bring almost zero values of M to zero to reduce dependencies
		snap_to_zero(M);

//print_lin_expressions(std::cout,M,M_dom);
		// Find which variables can be eliminated
		positional_vdomain wp_dom; // u+v variables that can be eliminated
		positional_vdomain wpp_dom; // u+v variables that remain
		i=0;
		unsigned int c = std::min(nb_cons,w_dom.size());
		unsigned int level = 0;
		while (i<w_dom.size() && level<w_dom.size() && level<M.size1()) {
			j=0;
			bool is_unit = true;
			while (j<M.size1() && is_unit) {
				if (j==level) {
					is_unit = is_unit && is_MEQ(M(j,i),scalar_type(1));
				} else {
					is_unit = is_unit && is_MEQ(M(j,i),scalar_type(0));
				}
				++j;
			}
			if (is_unit) {
				// add variable i to wp_dom
				wp_dom.add_variable(M_dom.get_variable(i));
				++level;
			} else {
				wpp_dom.add_variable(M_dom.get_variable(i));
			}
			++i;
		}
		// add the remaining variables to wpp
		while (i<w_dom.size()) {
			wpp_dom.add_variable(M_dom.get_variable(i));
			++i;
		}

		assert(w_dom.size()==wp_dom.size()+wpp_dom.size());

		// i is the number of u's and v's that can be eliminated
		size_t nb_w_elim = level;

		if (nb_w_elim>0) {
			LOGGER_OS(DEBUG6,__FUNCTION__) << "eliminated " << wp_dom.size() << "/" << w_dom.size() << " algebraic variables: " << wp_dom << std::endl;

			// retrieve the matrices
			size_t nb_w = w_dom.size();
			size_t nb_w_rem = wpp_dom.size();
			i=nb_w_elim;
			j=0;
			u=nb_w_rem;
			v=nb_w_elim;
//			std::cout << "assigning Dp" << std::endl << std::flush;
			matrix<scalar_type> Dp(nb_w_elim,wpp_dom.size());
			for (size_t q=0;q<wpp_dom.size();++q) {
				size_t M_index = M_dom.pos(wpp_dom.get_variable(q));
				for (size_t k=0;k<nb_w_elim;++k) {
					Dp(k,q)=M(k,M_index);
				}
			}
			i = nb_w;
			u = n;
			matrix<scalar_type> Cp = M.project_submatrix(j,j+v,i,i+u);
			matrix<scalar_type> Cpp = M.project_submatrix(j+v,nb_cons,i,i+u);
			i+=u;
			u = 1;
			matrix<scalar_type> dpm = M.project_submatrix(j,j+v,i,i+u);
			vector<scalar_type> dp = dpm.vector_from_column(0);
			matrix<scalar_type> dppm = M.project_submatrix(j+v,nb_cons,i,i+u);
			vector<scalar_type> dpp = dppm.vector_from_column(0);

			// construct B_w that applies to u and v
			matrix<scalar_type> B_w(n,nb_w);
			B_w.submatrix_assign(B,0,n,0,B.size2());

			matrix<scalar_type> Bp(n,wp_dom.size());
			for (i=0;i<wp_dom.size();++i) {
				size_t M_index = w_dom.pos(wp_dom.get_variable(i));
				Bp.assign_column(i,B_w.vector_from_column(M_index));
			}

			matrix<scalar_type> Bpp(n,wpp_dom.size());
			for (i=0;i<wpp_dom.size();++i) {
				size_t M_index = w_dom.pos(wpp_dom.get_variable(i));
				Bpp.assign_column(i,B_w.vector_from_column(M_index));
			}

//			std::cout << "elim: " << nb_w_elim << ", Bp: " << Bp << ", Cp: " << Cp << ", dp: " << dp << ", Dp " << Dp << std::endl << std::flush;
			A_elim = A - Bp*Cp;
			B_elim = Bpp - Bp*Dp; // over domain wpp_dom

			snap_to_zero(B_elim);

			LOGGER_OS(DEBUG6,__FUNCTION__) << "remaining " << wpp_dom.size() << "/" << w_dom.size() << " algebraic variables: " << wpp_dom << std::endl;
			//LOGGER_OS(DEBUG6,__FUNCTION__) << "input matrix: " << B_elim;

			alpha_elim = alpha - Bp*dp;
			// return the new domain for B
			B_dom = wpp_dom;

			// construct F_w that applies to u and v
			assert(slack_dom.size()==nb_ineq);
			matrix<scalar_type> F_w(slack_dom.size(),w_dom.size());
			// diagonal on the slack variables
			for (i=0;i<slack_dom.size();++i) {
				F_w(i,w_dom.pos(slack_dom.get_variable(i))) = scalar_type(1);
			}
			matrix<scalar_type> Fp(nb_ineq,wp_dom.size());
			for (i=0;i<wp_dom.size();++i) {
				size_t M_index = w_dom.pos(wp_dom.get_variable(i));
				Fp.assign_column(i,F_w.vector_from_column(M_index));
			}
			matrix<scalar_type> Fpp(nb_ineq,wpp_dom.size());
			for (i=0;i<wpp_dom.size();++i) {
				size_t M_index = w_dom.pos(wpp_dom.get_variable(i));
				vector<scalar_type> v = F_w.vector_from_column(M_index);
				Fpp.assign_column(i,v);
			}

			positional_vdomain F_dom = compose(x_dom,wpp_dom);
			matrix<scalar_type> red_F(nb_ineq,F_dom.size());
			matrix<scalar_type> Fx = Fp*Cp;
			matrix<scalar_type> Fwpp = Fp*Dp-Fpp;
			vector<scalar_type> f = Fp*dp;
			i = 0;
			u = Fx.size2();
			red_F.submatrix_assign(Fx,0,nb_ineq,i,i+u);
			i+=u;
			u=Fwpp.size2();
			red_F.submatrix_assign(Fwpp,0,nb_ineq,i,i+u);

//			std::cout << "final: A=" << A_elim << ", B=" << B_elim << ", b=" << alpha_elim << std::endl;
//			std::cout << "final: F=" << red_F << f << std::endl;


			if (!red_F.is_zero() || !f.is_zero()) {
				reduced_cons.push_back(construct_linear_constraint_system(red_F,-f,F_dom,ineq_signs));
//				std::cout << "F constraints: " << reduced_cons << std::endl;
			}
			if (!Cpp.is_zero() || !dpp.is_zero()) {
				// construct equality signs
				std::vector<sign> eq_signs(Cpp.size1(),EQ);
				reduced_cons.push_back(construct_linear_constraint_system(Cpp,-dpp,x_dom,eq_signs));
//				std::cout << "F+eq constraints: " << reduced_cons << std::endl;
			}
			if (!reduced_cons.empty()) {
				remove_unused_variables(reduced_cons);
//				std::cout << "w/o unused: " << reduced_cons << std::endl;
				reduced_cons.collapse_inequalities();
//				std::cout << "collapsed: " << reduced_cons << std::endl;
				LOGGER_OS(DEBUG6,__FUNCTION__) << "reduced domain from " << F_dom.size() << " to " << reduced_cons.get_variable_ids().size() << std::endl;
			}


//			std::cout << "constraints: " << reduced_cons << std::endl;
		} else {
			throw std::runtime_error(
					"Alg_var_elim: Elimination: Could not solve algebraic equations");
		}
	}
	// For now we don't do any processing of inequality constraints, so let's just hand these back
	reduced_cons.push_back(disregarded_cons);
}
template<typename scalar_type>
void alg_var_elim<scalar_type>::get_cons_matrices(matrix<scalar_type>&C_mat, matrix<scalar_type>&D_mat,matrix<scalar_type>&E_mat) const{
	C_mat = C;
	D_mat = D;
	E_mat = E;
}
template<typename scalar_type>
void alg_var_elim<scalar_type>::print() const{
	std::cout << "C:\n" << C << std::endl;
	std::cout << "D:\n" << D << std::endl;
	std::cout << "E:\n" << E << std::endl;
	std::cout << "beta:" << beta << std::endl;
	std::cout << "disregarded constraints:" << disregarded_cons << std::endl;
}

}
