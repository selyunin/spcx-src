/*
 * omega_interpol.h
 *
 *  Created on: Jan 20, 2011
 *      Author: frehse
 */

#ifndef OMEGA_INTERPOL_H_
#define OMEGA_INTERPOL_H_

#include "omega_model.h"
#include "core/continuous/polyhedra/hyperbox/finite_hyperbox_utility.h"

namespace continuous {
namespace support_function {

/** Overapproximation of the reachable set in a time interval,
 * using a forward error model. */

template<typename scalar>
class omega_interpol: public omega_model<scalar> {
public:
	typedef omega_model<scalar> base_class;
	typedef typename base_class::model_ptr model_ptr;

	typedef math::vdom_matrix<scalar> matrix_type;
	typedef math::vdom_vector<scalar> vector_type;

	/** Construct the model with given sets.
	 *
	 * X0 and U are passed as pointers so they can be null.
	 */
	omega_interpol(const matrix_type& A, const scalar& delta,
			const support_function_provider* X0,
			const support_function_provider* U);

	/** Construct the model with a state_handler to reduce redundant computations.
	 */
	omega_interpol(model_ptr state_handler);

	/** Create a full clone.
	 * */
	omega_interpol<scalar>* clone();

	/** Create a clone that shares all state except the parts that depend on l.
	 * */
	omega_interpol<scalar>* client();

	using base_class::get_l;
	using base_class::get_l_next;

protected:
	/** Update delta */
	void update_impl();

	/** Obtain e^{A\delta}^T for current delta, implementation */
	const matrix_type& exp_AdeltaT_impl();

	/** Compute the support of omega w.r.t. direction l and the support value of X0 and U in direction l.
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

	/** Compute the support of omega w.r.t. direction l and the support value of X0 and U in direction l.
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
	omega_interpol<scalar>* handler();

	/** Create a map from a matrix */
	math::affine_map<scalar> to_map(const matrix_type& M);

	/** Compute error model for given A,delta,X0,U. */
	void compute_err_model();

	// Parameters that change with delta
	matrix_type Phi; //=e^(A*delta)
	matrix_type Phi_T; //=Phi^T
	matrix_type Phi1; //=A^{-1}(e^(A*delta)-I)
	matrix_type Phi2; //=A^{-2}(e^(A*delta)-I-delta*A)

	typedef finite_hyperbox<scalar> box_type;
	box_type E_X0;
	box_type E_U;
};

template<typename scalar>
omega_interpol<scalar>* omega_interpol<scalar>::handler() {
	return static_cast<omega_interpol<scalar>*> (omega_model<scalar>::handler());
}

template<typename scalar>
omega_interpol<scalar>* omega_interpol<scalar>::clone() {
	return new omega_interpol<scalar> (*this);
}

template<typename scalar>
omega_interpol<scalar>* omega_interpol<scalar>::client() {
	return new omega_interpol<scalar> (get_ptr());
}

template<typename scalar>
omega_interpol<scalar>::omega_interpol(const matrix_type& A,
		const scalar& delta, const support_function_provider* X0,
		const support_function_provider* U) :
	omega_model<scalar> (A, delta, X0, U) {
	update_impl();
}

template<typename scalar>
omega_interpol<scalar>::omega_interpol(model_ptr state_handler) :
	omega_model<scalar> (state_handler) {
}

template<typename scalar>
scalar omega_interpol<scalar>::Omega_support_impl(scalar& err) {
	// This model is taken from Colas' thesis, p. 77/78
	scalar a(0), b(0), c(0), d(0), e(0);
	scalar rho_omega;
	if (get_U()) {
		d = get_delta() * get_rho_U();
		e = rho(get_l(), handler()->E_U);
	}
	if (get_X0()) {
		a = get_rho_X0();
		b = get_rho_X0_next();
		c = rho(get_l(), handler()->E_X0);
		// get the max of the quadratic function on the interval 0<=lambda<=1
		scalar lambda_max;
		rho_omega = quadratic_max<scalar> (e - c, -a + b + c + d, a, scalar(0),
				scalar(1), lambda_max);
		err = lambda_max * (scalar(1) - lambda_max) * c + lambda_max
				* lambda_max * e;
	} else {
		err = e;
		rho_omega = d + e;
	}
	return rho_omega;
}

template<typename scalar>
scalar omega_interpol<scalar>::Psi_support_impl(scalar& err) {
	if (get_U()) {
		err = rho(get_l(), handler()->E_U);
		return get_delta() * get_rho_U() + err;
	} else {
		err = scalar(0);
		return scalar(0);
	}
}

template<typename scalar>
void omega_interpol<scalar>::update_impl() {
	/* New computation of Phi at the same time as
	 * special matrices
	 * Phi=e^(A*delta)
	 * Phi1=A^{-1}(e^(A*delta)-I)
	 * Phi2=A^{-2}(e^(A*delta)-I-delta*A)
	 */
	positional_vdomain dom = get_A().domain();
	Phi = math::matrix_exponential(get_delta()*get_A());
	Phi_T = Phi.transpose(); // Phi_T is the transpose of Phi.

	matrix_type A_abs = get_A().abs();
	matrix_type Phi_abs = matrix_type(dom, dom); // not used
	matrix_type Phi1_abs = matrix_type(dom, dom);
	matrix_type Phi2_abs = matrix_type(dom, dom);
	matrix_type M2; // = A-A*exp(\delta*A)
	matrix_type M4; // = A^2*exp(\delta*A)
	matrix_type Mtemp;

	math::get_special_matrices(A_abs.get_matrix(), get_delta(),
			Phi_abs.get_matrix(), Phi1_abs.get_matrix(), Phi2_abs.get_matrix());

	Mtemp = get_A() * Phi;
	M2 = get_A() - Mtemp;
	M4 = get_A() * Mtemp;

	if (get_X0()) {
		box_type S1 = finite_symmetric_bounding_box(*get_X0(), to_map(M2));
		box_type BX1 = finite_bounding_box(S1, to_map(Phi1_abs));
		box_type S2 = finite_symmetric_bounding_box(*get_X0(), to_map(M4));
		box_type BX2 = finite_bounding_box(S2, to_map(Phi2_abs));
		E_X0 = BX1 + BX2;
	} else {
		vector_type zeros(dom);
		E_X0 = finite_hyperbox<scalar> (zeros, zeros);
	}

	if (get_U()) {
		box_type SU = finite_symmetric_bounding_box(*get_U(), to_map(get_A()));
		E_U = finite_bounding_box(SU, to_map(Phi2_abs));
	} else {
		vector_type zeros(dom);
		E_U = finite_hyperbox<scalar> (zeros, zeros);
	}

	//std::cout << S1 << BX1 << std::endl << S2 << BX2 << std::endl << E_X0 << std::endl;
}

template<typename scalar>
const typename omega_interpol<scalar>::matrix_type& omega_interpol<scalar>::exp_AdeltaT_impl() {
	return Phi_T;
}

template<class T>
math::affine_map<T> omega_interpol<T>::to_map(const matrix_type& M) {
	//	positional_vdomain dom = get_A().domain();
	//	typedef typename math::affine_map<T>::vdom_matrix_type dom_matrix;
	//	return math::affine_map<T>(dom_matrix(dom, dom, M));
	return math::affine_map<T>(M);
}

}
}

#endif /* OMEGA_INTERPOL_H_ */
