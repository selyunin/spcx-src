/*
 * omega_forward.h
 *
 *  Created on: Jan 19, 2011
 *      Author: frehse
 */

#ifndef OMEGA_FORWARD_H_
#define OMEGA_FORWARD_H_

#include "omega_model.h"

namespace continuous {
namespace support_function {

/** Overapproximation of the reachable set in a time interval,
 * using a forward error model. */

template<typename scalar>
class omega_forward: public omega_model<scalar> {
public:
	typedef omega_model<scalar> base_class;
	typedef typename base_class::model_ptr model_ptr;

	typedef math::vdom_matrix<scalar> matrix_type;
	typedef math::vdom_vector<scalar> vector_type;

	/** Construct the model with given sets.
	 *
	 * X0 and U are passed as pointers so they can be null.
	 */
	omega_forward(const matrix_type& A, const scalar& delta,
			const support_function_provider* X0,
			const support_function_provider* U);

	/** Construct the model with a state_handler to reduce redundant computations.
	 */
	omega_forward(model_ptr state_handler);

	/** Create a full clone.
	 * */
	omega_forward<scalar>* clone();

	/** Create a clone that shares all state except the parts that depend on l.
	 * */
	omega_forward<scalar>* client();

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
	omega_forward<scalar>* handler();

	/**
	 * beta = mu * (e^{delta*||A||} - 1 - delta*||A||)/ ||A||
	 */
	void compute_beta();

	/**
	 * alpha = (e^{delta * ||A||} - 1 - delta*||A||) * gamma
	 */
	void compute_alpha();

	/**
	 * gamma = ||X0|| + ||U||/||A||
	 */
	void compute_gamma();

	scalar mu;
	scalar norm_A;
	scalar gamma;

	// Parameters that change with delta
	matrix_type Phi; //=e^(A*delta)
	matrix_type Phi_T; //=Phi^T
	matrix_type Phi1; //=A^{-1}(e^(A*delta)-I)
	matrix_type Phi2; //=A^{-2}(e^(A*delta)-I-delta*A)
	scalar alpha;
	scalar beta;

	/** Computes the support function on a ball of radius r.
	 *
	 * Computes mx l.x, x \in {|x| <= r},
	 * where |x| is the infinity norm.
	 * \param l vector
	 * \param r radius of the ball
	 */
	scalar ball_sf(const scalar& r, const math::vector<scalar>& l);
};

template<typename scalar>
omega_forward<scalar>* omega_forward<scalar>::handler() {
	return static_cast<omega_forward<scalar>*> (omega_model<scalar>::handler());
}

template<typename scalar>
omega_forward<scalar>* omega_forward<scalar>::clone() {
	return new omega_forward<scalar> (*this);
}

template<typename scalar>
omega_forward<scalar>* omega_forward<scalar>::client() {
	return new omega_forward<scalar> (get_ptr());
}

template<class T>
T omega_forward<T>::ball_sf(const T& r, const math::vector<T>& l) {
	/* the support function of a ball of radius r in direction l
	 * is r ||l||_1 if it is a ball for the infinity norm
	 * and r ||l||_\infty if it is a ball for the 1-norm. */
	T u = l.one_norm(); // norm should return one norm, since matrix is in infinity-norm
	return u * r;
}

template<typename scalar>
omega_forward<scalar>::omega_forward(const matrix_type& A, const scalar& delta,
		const support_function_provider* X0, const support_function_provider* U) :
	omega_model<scalar> (A, delta, X0, U) {
	norm_A = get_A().infinity_norm();
	if (get_U())
		get_max_infinity_norm(*get_U(), mu);

	compute_gamma();
	update_impl();
}

template<typename scalar>
omega_forward<scalar>::omega_forward(model_ptr state_handler) :
	omega_model<scalar> (state_handler) {
}

template<typename scalar>
scalar omega_forward<scalar>::Omega_support_impl(scalar& err) {
	/* Compute the support function of l on chull(X_0,\Phi X_0 + delta U + alpha Ball).
	 * Use \rho_l(chull(P,Q))=max(\rho_l(P),\rho_l(Q)). */
	err = ball_sf(handler()->alpha, get_l().get_vector());
	scalar U_term(0);
	if (get_U())
		U_term = get_delta() * get_rho_U();
	return std::max(get_rho_X0(), get_rho_X0_next() + U_term + err);
}

template<typename scalar>
scalar omega_forward<scalar>::Psi_support_impl(scalar& err) {
	err = ball_sf(handler()->beta, get_l().get_vector());
	if (get_U())
		return get_delta() * get_rho_U() + err;
	else
		return scalar(0);
}

template<typename scalar>
void omega_forward<scalar>::update_impl() {
	compute_beta();
	compute_alpha();

	/* New computation of Phi at the same time as
	 * special matrices
	 * Phi=e^(A*delta)
	 * Phi1=A^{-1}(e^(A*delta)-I)
	 * Phi2=A^{-2}(e^(A*delta)-I-delta*A)
	 */
	positional_vdomain dom = get_A().domain();
	Phi = matrix_type(dom, dom);
	Phi1 = matrix_type(dom, dom);
	Phi2 = matrix_type(dom, dom);
	math::get_special_matrices(get_A().get_matrix(), get_delta(),
			Phi.get_matrix(), Phi1.get_matrix(), Phi2.get_matrix());
	Phi_T = Phi.transpose(); // Phi_T is the transpose of Phi.
}

template<typename scalar>
const typename omega_forward<scalar>::matrix_type& omega_forward<scalar>::exp_AdeltaT_impl() {
	return Phi_T;
}

template<class T> void omega_forward<T>::compute_beta() {
	if (norm_A == T(0)) {
		beta = T(0);
	} else {
		double a, b;
		b = convert_element<double> (norm_A * get_delta());
		a = std::pow(M_E, b);

		beta = (T(a) - T(1) - get_delta() * norm_A) * mu / norm_A;
	}
}

template<class T> void omega_forward<T>::compute_alpha() {
	double a, b;
	b = convert_element<double> (norm_A * get_delta());
	a = std::pow(M_E, b);
	alpha = (T(a) - T(1) - get_delta() * norm_A) * gamma;
}

template<class T> void omega_forward<T>::compute_gamma() {
	if (get_X0())
		get_max_infinity_norm(*get_X0(), gamma);
	else
		gamma = T(0);
	if (norm_A > T(0)) {
		gamma += mu / norm_A;
	} else {
		// if A_norm is zero we get a NaN otherwise
		gamma = T(0);
	}
}

}
}

#endif /* OMEGA_FORWARD_H_ */
