/*
 * omega_interpol_fb_fb.h
 *
 *  Created on: Jan 20, 2011
 *      Author: frehse
 */

#ifndef OMEGA_INTERPOL_FB_H_
#define OMEGA_INTERPOL_FB_H_

#include "omega_model.h"
#include "core/continuous/polyhedra/hyperbox/finite_hyperbox_utility.h"
#include "core/continuous/support_function/sf_base/sf_unary_ref.h"

namespace continuous {
namespace support_function {

/** Overapproximation of the reachable set in a time interval,
 * using a forward error model.
 *
 * By default, the error bounds are
 * e = max_{\lambda in [0,1]} \rho(Omega_\lambda) - \rho(Reach_\lambda).
 *
 * To use the error of the actual computed constraint,
 * set use_error_at_t = false;
 * */

template<typename scalar>
class omega_interpol_fb: public omega_model<scalar> {
public:
	typedef omega_model<scalar> base_class;
	typedef typename base_class::model_ptr model_ptr;

	typedef math::vdom_matrix<scalar> matrix_type;
	typedef math::vdom_vector<scalar> vector_type;

	/** Construct the model with given sets.
	 *
	 * X0 and U are passed as pointers so they can be null.
	 * If minerror is set to true, the error model return
	 * the difference between the largest and the smallest
	 * support function value of Omega (conservative
	 * estimate that accounts for approx. error).
	 */
	omega_interpol_fb(const matrix_type& A, const scalar& delta,
			const support_function_provider* X0,
			const support_function_provider* U, bool minerror = false);

	/** Construct the model with a state_handler to reduce redundant computations.
	 *
	 * If the state handler is of type omega_interpol_fb, all parameters
	 * (minerror) are copied.
	 */
	omega_interpol_fb(model_ptr state_handler);

	/** Create a full clone.
	 * */
	omega_interpol_fb<scalar>* clone();

	/** Create a clone that shares all state except the parts that depend on l.
	 * */
	omega_interpol_fb<scalar>* client();

	using base_class::get_l;
	using base_class::get_l_next;

protected:
	/** Update delta */
	void update_impl();

	/** Obtain e^{A\delta}^T for current delta, implementation */
	const matrix_type& exp_AdeltaT_impl();

	/** Compute the support of Omega w.r.t. direction l.
	 *
	 * The overapproximation error is returned in err.
	 *
	 * The implentation can use the following member variables:
	 * - my_l           : current l
	 * - my_lnext       : l^Te^{A\delta}
	 * - my_rho_X0      : rho(l,X0)
	 * - my_rho_X0_next : rho(l,e^{At}X0)
	 * - my_rho_U       : rho(l,U)
	 *
	 * These values are defined if X0, respectively U are not null.
	 * The implementation must take null X0 and U into account. */
	scalar Omega_support_impl(scalar& err);

	/** Compute the support of Psi w.r.t. direction l.
	 *
	 * The overapproximation error is returned in err.
	 *
	 * The implentation can use the following member variables:
	 * - my_l           : current l
	 * - my_lnext       : l^Te^{A\delta}
	 * - my_rho_X0      : rho(l,X0)
	 * - my_rho_X0_next : rho(l,e^{At}X0)
	 * - my_rho_U       : rho(l,U)
	 *
	 * These values are defined if X0, respectively U are not null.
	 * The implementation must take null X0 and U into account. */
	scalar Psi_support_impl(scalar& err);

	/** Compute the support of Omega w.r.t. direction l at lambda by interpolation.
	 *
	 * Returns the support of Omega that covers Reach_[t+lambda1*delta,t+lambda2*delta].
	 */
	scalar Omega_support_interp_impl(const scalar& lambda1,
			const scalar& lambda2, scalar& err);

	/** Compute the support of Psi w.r.t. direction l at lambda by interpolation.
	 *
	 * Returns the support of Psi that covers Reach_[t+lambda*delta,t+lambda*delta].
	 */
	scalar Psi_support_interp_impl(const scalar& lambda, scalar& err);

	using base_class::get_ptr;
	using base_class::get_A;
	using base_class::get_delta;
	using base_class::get_X0;
	using base_class::get_U;
	using base_class::get_rho_X0;
	using base_class::get_rho_U;
	using base_class::get_rho_X0_next;

private:
	/** Get a pointer to the handler */
	omega_interpol_fb<scalar>* handler();

	/** Create a map from a matrix */
	math::affine_map<scalar> to_map(const matrix_type& M);

	/** Compute error model for given A,delta,X0,U. */
	void compute_err_model();

	/** Optimized equation */
	scalar product_lambda_mins(const vector_type& ef, const vector_type& eb, const vector_type& l,const scalar& lambda1);

	// Parameters that change with delta
	matrix_type Phi; //=e^(A*delta)
	matrix_type Phi_T; //=Phi^T
	matrix_type Phi1; //=A^{-1}(e^(A*delta)-I)
	matrix_type Phi2; //=A^{-2}(e^(A*delta)-I-delta*A)

	typedef finite_hyperbox<scalar> box_type;
	box_type E_X0f;
	box_type E_X0b;
	box_type E_U;
	support_function_provider::ptr negAPhi2U;

	bool my_minerror; // if true, use difference between max and min Omega_t as error
};

template<typename scalar>
omega_interpol_fb<scalar>* omega_interpol_fb<scalar>::handler() {
	return static_cast<omega_interpol_fb<scalar>*> (omega_model<scalar>::handler());
}

template<typename scalar>
omega_interpol_fb<scalar>* omega_interpol_fb<scalar>::clone() {
	return new omega_interpol_fb<scalar> (*this);
}

template<typename scalar>
omega_interpol_fb<scalar>* omega_interpol_fb<scalar>::client() {
	return new omega_interpol_fb<scalar> (get_ptr());
}

template<typename scalar>
omega_interpol_fb<scalar>::omega_interpol_fb(const matrix_type& A,
		const scalar& delta, const support_function_provider* X0,
		const support_function_provider* U, bool minerror) :
	omega_model<scalar> (A, delta, X0, U), my_minerror(minerror) {
	update_impl();
}

template<typename scalar>
omega_interpol_fb<scalar>::omega_interpol_fb(model_ptr state_handler) :
	omega_model<scalar> (state_handler) {
	// copy parameters if it's an omega_interpol_fb
	if (omega_interpol_fb* p =  dynamic_cast<omega_interpol_fb*>(state_handler.get())) {
		my_minerror = p->my_minerror;
	}
}

template<typename scalar>
scalar omega_interpol_fb<scalar>::Omega_support_impl(scalar& err) {
	return Omega_support_interp_impl(scalar(0), scalar(1), err);
}

template<typename scalar>
scalar omega_interpol_fb<scalar>::Psi_support_impl(scalar& err) {
	return Psi_support_interp_impl(scalar(1), err);
}

template<typename scalar>
scalar omega_interpol_fb<scalar>::product_lambda_mins(const vector_type& ef, const vector_type& eb, const vector_type& l,const scalar& lambda1) {
	// computes scalar_product(math::min(lambda1 * ef, (scalar(1) - lambda1) * eb),  math::vec_abs(l));
	scalar x = scalar(0);
	scalar one_m_lambda1 = scalar(1) - lambda1;
	for (size_t i = 0; i < ef.size(); ++i) {
		scalar lambda1_ef = abs(lambda1 * ef[i]);
		scalar one_lambda1_eb = abs(one_m_lambda1 * eb[i]);
		x += std::min(lambda1_ef,one_lambda1_eb) * abs(l[i]);
	}
	return x;
}

template<typename scalar>
scalar omega_interpol_fb<scalar>::Omega_support_interp_impl(
		const scalar& lambda1, const scalar& lambda2, scalar& err) {
	// This model is taken from Colas' thesis, p. 77/78
	scalar a(0), b(0), c(0), d(0), e(0);
	scalar ufactor(0); // rho_negAPhi2U
	scalar rho_omega;
	scalar cf(0), cb(0);
	scalar rho_omegaf, errf; // forward
	scalar rho_omegab, errb; // backward
	scalar rho_omegai, erri; // interpolation
	scalar rho_omegamin; // for min error computation
	if (get_U()) {
		d = get_delta() * get_rho_U();
		e = rho(get_l(), handler()->E_U);
		ufactor = rho(get_l(), *handler()->negAPhi2U);
	}
	if (get_X0()) {
		a = get_rho_X0();
		b = get_rho_X0_next();

		/** This is the proper version from the CAV11 paper */
		//vector_type l_abs = math::vec_abs(get_l());
		vector_type ef = math::vec_abs(rho_vec(get_l(), handler()->E_X0f));
		vector_type eb = math::vec_abs(rho_vec(get_l(), handler()->E_X0b));

		// lambda = lambda1
//		cf = scalar_product(math::min(lambda1 * ef, (scalar(1) - lambda1) * eb), l_abs);
		cf = product_lambda_mins(ef,eb,get_l(),lambda1);
		rho_omega = e * lambda1 * lambda1 + (-a + b + d) * lambda1 + cf + a;
		//std::cout << "at lambda="<< lambda_switch << " -> " << rho_omegaf << ", cf=" << cf << " min=" << math::min(lambda_switch*ef,(scalar(1)-lambda_switch)*eb) << " labs="<<l_abs <<std::endl;
		err = lambda1 * (ufactor + cf) + e * lambda1 * lambda1;
		rho_omegamin = rho_omega - err;

		// lambda = lambda2
		//cf = scalar_product(math::min(lambda2 * ef, (scalar(1) - lambda2) * eb), l_abs);
		cf = product_lambda_mins(ef,eb,get_l(),lambda2);
		rho_omegaf = e * lambda2 * lambda2 + (-a + b + d) * lambda2 + cf + a;
		//std::cout << "at lambda="<< lambda_switch << " -> " << rho_omegaf << ", cf=" << cf << " min=" << math::min(lambda_switch*ef,(scalar(1)-lambda_switch)*eb) << " labs="<<l_abs <<std::endl;
		errf = lambda2 * (ufactor + cf) + e * lambda2 * lambda2;
		if (rho_omegaf > rho_omega) {
			rho_omega = rho_omegaf;
		}
		if (errf > err) {
			err = errf;
		}
		if (rho_omegaf - errf < rho_omegamin)
			rho_omegamin = rho_omegaf - errf;

		// cycle through switching lambdas
		for (size_t i = 0; i < ef.size(); ++i) {
			if (math::definitely(math::numeric::is_GT(ef[i] + eb[i], scalar(0)))) {
				scalar lambda_switch = eb[i] / (ef[i] + eb[i]);
				if (lambda_switch > lambda1 && lambda_switch < lambda2) {
					//cf = scalar_product(math::min(lambda_switch * ef,(scalar(1) - lambda_switch) * eb), l_abs);
					cf = product_lambda_mins(ef,eb,get_l(),lambda_switch);
					rho_omegaf = e * lambda_switch * lambda_switch + (-a + b
							+ d) * lambda_switch + cf + a;
					//std::cout << "at lambda="<< lambda_switch << " -> " << rho_omegaf << ", cf=" << cf << " min=" << math::min(lambda_switch*ef,(scalar(1)-lambda_switch)*eb) << " labs="<<l_abs <<std::endl;
					errf = lambda_switch * (ufactor + cf) + e * lambda_switch
							* lambda_switch;
					if (rho_omegaf > rho_omega) {
						rho_omega = rho_omegaf;
					}
					if (errf > err) {
						err = errf;
					}
					if (rho_omegaf - errf < rho_omegamin)
						rho_omegamin = rho_omegaf - errf;
				}
			}
		}
	} else {
		// lambda = lambda1
		rho_omega = e * lambda1 * lambda1 + d * lambda1;
		err = lambda1 * (ufactor) + e * lambda1 * lambda1;
		rho_omegamin = rho_omega - err;

		// lambda = lambda2
		rho_omegaf = e * lambda2 * lambda2 + d * lambda2;
		errf = lambda2 * (ufactor) + e * lambda2 * lambda2;
		if (rho_omegaf > rho_omega) {
			rho_omega = rho_omegaf;
		}
		if (errf > err) {
			err = errf;
		}
		if (rho_omegaf - errf < rho_omegamin)
			rho_omegamin = rho_omegaf - errf;
	}
	if (my_minerror)
		err = rho_omega - rho_omegamin;
	return rho_omega;
}

template<typename scalar>
scalar omega_interpol_fb<scalar>::Psi_support_interp_impl(const scalar& lambda,
		scalar& err) {
	if (get_U()) {
		err = scalar(0);
		scalar err_psi = rho(get_l(), handler()->E_U);
		scalar rho_psi = lambda * get_delta() * get_rho_U() + err_psi;

		// compute error
		// err = err_psi + rho_(-A Phi_2 U)
		scalar ufactor = rho(get_l(), *handler()->negAPhi2U);
		err = err_psi + ufactor;

		return rho_psi;
	} else {
		err = scalar(0);
		return scalar(0);
	}
}

template<typename scalar>
void omega_interpol_fb<scalar>::update_impl() {
	/* New computation of Phi at the same time as
	 * special matrices
	 * Phi=e^(A*delta)
	 * Phi1=A^{-1}(e^(A*delta)-I)
	 * Phi2=A^{-2}(e^(A*delta)-I-delta*A)
	 */
	const matrix_type& A = get_A();
	positional_vdomain dom = A.domain();
	Phi = math::matrix_exponential(get_delta() * get_A());

	Phi_T = Phi.transpose(); // Phi_T is the transpose of Phi.

	matrix_type A_abs = A.abs();
	matrix_type Phi_abs = matrix_type(dom, dom); // not used
	matrix_type Phi1_abs = matrix_type(dom, dom);
	matrix_type Phi2_abs = matrix_type(dom, dom);
	matrix_type M2; // = A-A*exp(\delta*A)
	matrix_type M4; // = A^2*exp(\delta*A)
	matrix_type Mtemp;

	math::get_special_matrices(A_abs.get_matrix(), get_delta(),
			Phi_abs.get_matrix(), Phi1_abs.get_matrix(), Phi2_abs.get_matrix());

	Mtemp = A * Phi;
	M2 = A - Mtemp;
	M4 = A * Mtemp;

	if (get_X0()) {
		//		box_type S1 = finite_symmetric_bounding_box(*get_X0(), to_map(M2));
		//		box_type BX1 = finite_bounding_box(S1, to_map(Phi1_abs));
		//		box_type S2 = finite_symmetric_bounding_box(*get_X0(), to_map(M4));
		//		box_type BX2 = finite_bounding_box(S2, to_map(Phi2_abs));
		//		E_X0 = BX1 + BX2;
		continuous::support_function_provider const& X0 = *get_X0();
		math::affine_map<scalar> A2_map = to_map(A * A);
		box_type S2f = finite_symmetric_bounding_box(X0, A2_map);
		//std::cout << "Phi2_abs:" << Phi2_abs << std::endl;
		box_type BX2f = finite_bounding_box(S2f, to_map(Phi2_abs));
		//std::cout << "S2f:" << S2f << std::endl;
		//std::cout << "BX2f:" << BX2f << std::endl;
		box_type S2b = finite_symmetric_bounding_box(X0, to_map(M4));
		box_type BX2b = finite_bounding_box(S2b, to_map(Phi2_abs));
		//std::cout << "S2b:" << S2b << std::endl;
		//std::cout << "BX2b:" << BX2b << std::endl;
		E_X0f = BX2f;
		E_X0b = BX2b;
	} else {
		vector_type zeros(dom);
		E_X0f = finite_hyperbox<scalar> (zeros, zeros);
		E_X0b = finite_hyperbox<scalar> (zeros, zeros);
	}

	if (get_U()) {
		box_type SU = finite_symmetric_bounding_box(*get_U(), to_map(get_A()));
		E_U = finite_bounding_box(SU, to_map(Phi2_abs));

		// we need Phi2 for the error bounds
		matrix_type Phi1 = matrix_type(dom, dom); // not used
		matrix_type Phi2 = matrix_type(dom, dom);
		math::get_special_matrices(A.get_matrix(), get_delta(),
				Phi.get_matrix(), Phi1.get_matrix(), Phi2.get_matrix());
		negAPhi2U = support_function_provider::ptr(new sf_unary_ref<scalar> (
				*get_U(), to_map(-A * Phi2)));
	} else {
		vector_type zeros(dom);
		E_U = finite_hyperbox<scalar> (zeros, zeros);
		negAPhi2U = support_function_provider::ptr();
	}

	//	std::cout << "U:" << handler()->E_U <<  handler()->E_U.is_empty() << std::endl;
	//	std::cout << "X0f:" << handler()->E_X0f <<  handler()->E_X0f.is_empty() << std::endl;
	//	std::cout << "X0b:" << handler()->E_X0b <<  handler()->E_X0b.is_empty() << std::endl << std::flush;
	//std::cout << S1 << BX1 << std::endl << S2 << BX2 << std::endl << E_X0 << std::endl;
}

template<typename scalar>
const typename omega_interpol_fb<scalar>::matrix_type& omega_interpol_fb<scalar>::exp_AdeltaT_impl() {
	return Phi_T;
}

template<class T>
math::affine_map<T> omega_interpol_fb<T>::to_map(const matrix_type& M) {
	//	positional_vdomain dom = get_A().domain();
	//	typedef typename math::affine_map<T>::vdom_matrix_type dom_matrix;
	//	return math::affine_map<T>(dom_matrix(dom, dom, M));
	return math::affine_map<T>(M);
}

}
}

#endif /* OMEGA_INTERPOL_FB_H_ */
