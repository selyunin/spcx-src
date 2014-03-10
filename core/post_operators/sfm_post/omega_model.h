/*
 * omega_model.h
 *
 *  Created on: Jan 19, 2011
 *      Author: frehse
 */

#ifndef OMEGA_MODEL_H_
#define OMEGA_MODEL_H_

#include "math/matrix.h"
#include "math/vector.h"
#include "math/vdom/vdom_matrix.h"
#include "math/numeric/comp.h"
#include "math/unique_scalar_to_value_store.h"
#include "core/continuous/support_function_provider.h"

namespace continuous {
namespace support_function {

/** Overapproximation of the reachable set in a time interval. */

template<typename scalar>
class omega_model: public boost::enable_shared_from_this<omega_model<scalar> > {
public:
	typedef math::vdom_matrix<scalar> matrix_type;
	typedef math::vdom_vector<scalar> vector_type;
	typedef boost::shared_ptr<omega_model<scalar> > model_ptr;
	typedef math::unique_scalar_to_value_store<scalar, model_ptr>
			model_cache_type;
	typedef boost::shared_ptr<model_cache_type> model_cache_ptr;

	/** Virtual destructor. */
	virtual ~omega_model();

	/** Create a full clone.
	 * */
	virtual omega_model<scalar>* clone() = 0;

	/** Create a clone that shares all state except the parts that depend on l.
	 * */
	virtual omega_model<scalar>* client() = 0;

	/** Return a shared_ptr to *this. */
	model_ptr get_ptr();

	/** Compute the support of X0 w.r.t. direction l. */
	virtual scalar rho_X0(const vector_type& l);

	/** Compute the support of U w.r.t. direction l. */
	virtual scalar rho_U(const vector_type& l);

	/** Compute the support of Omega w.r.t. direction l.
	 *
	 * The overapproximation error is returned in err. */
	virtual scalar Omega_support(const vector_type& l, scalar& err);

	/** Compute the support of Psi w.r.t. direction l.
	 *
	 * The overapproximation error is returned in err. */
	virtual scalar Psi_support(const vector_type& l, scalar& err);

	/** Initialize a sequence of computations of the support w.r.t. direction l.
	 *
	 * The sequence is computed with support_sequence_Omega
	 * and support_sequence_Psi, then support_sequence_next
	 * to advance to the next element. The model internally
	 * keeps track of the time starting with t0. */
	virtual void support_sequence_ini(const vector_type& l, const scalar& t0 =
			scalar(0));

	/** Compute the support of Omega in the sequence.
	 *
	 * The overapproximation error is returned in err. */
	virtual scalar support_sequence_Omega(scalar& err);

	/** Compute the support of Omega in the sequence.
	 *
	 * The overapproximation error is returned in err. */
	virtual scalar support_sequence_Psi(scalar& err);

	/** Advance to the next element of the sequence.
	 *
	 * The sequence needs to be initialized with support_sequence_ini first. */
	virtual void support_sequence_next();

	/** Returns whether the model can interpolate
	 *
	 * To be overloaded by derived classes that can interpolate
	 * The default implementation is to return false. */
	virtual bool support_sequence_interpolates() const;

	/** Compute the support of Omega between t1 and t2 by interpolation.
	 *
	 * The overapproximation error is returned in err. */
	virtual scalar support_sequence_Omega_support_interp(const scalar& t1,
			const scalar& t2, scalar& err);

	/** Compute the support of Psi at t by interpolation.
	 *
	 * The overapproximation error is returned in err. */
	virtual scalar support_sequence_Psi_support_interp(const scalar& t,
			scalar& err);

	/** Update delta */
	void update(const scalar& delta);

	/** Obtain e^{A\delta}^T for current delta */
	const matrix_type& exp_AdeltaT();

	/** Obtain current l in a sequence. */
	virtual const vector_type& get_l() const;

	/** Obtain e^{A\delta}^Tl for current delta and l in a sequence. */
	virtual const vector_type& get_l_next() const;

	/** Get the current value of delta. */
	const scalar& get_delta() const;

	/** Get the current value of t in a sequence. */
	const scalar& get_t() const;

	/** Get the next value of t in a sequence. */
	scalar get_t_next() const;

protected:
	/** Construct the model with given sets.
	 *
	 * X0 and U are passed as pointers so they can be null.
	 * If state_handler is not null, all updates will
	 * be handled by the state_handler to reduce redundant computations.
	 */
	omega_model(const matrix_type& A, const scalar& delta,
			const support_function_provider* X0,
			const support_function_provider* U);

	/** Construct the model with a state_handler to reduce redundant computations.
	 *
	 * This is a shallow copy and can be used to create models that have different
	 * l but the share same A, X0, U.
	 *
	 * Uses a pointer instead of reference so the signature doesn't override
	 * the copy constructor.
	 */
	omega_model(model_ptr state_handler);

	/** Update delta in derived classes */
	virtual void update_impl() = 0;

	/** Obtain e^{A\delta}^T for current delta, implementation */
	virtual const matrix_type& exp_AdeltaT_impl() = 0;

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
	virtual scalar Omega_support_impl(scalar& err) = 0;

	/** Compute the support of omega w.r.t. direction l and the support value of X0 and U in direction l.
	 *
	 * The overapproximation error is returned in err.
	 *
	 * The implentation can use the following member variables:
	 * - my_l           : current l
	 * - my_lnext       : l^Te^{A\delta}
	 * - my_rho_X0      :rho(l,X0)
	 * - my_rho_X0_next :rho(l,e^{At}X0)
	 * - my_rho_U       :rho(l,U)
	 *
	 * These values are defined if X0, respectively U are not null.
	 * The implementation must take null X0 and U into account. */
	virtual scalar Psi_support_impl(scalar& err) = 0;

	/** Compute the support of Omega w.r.t. direction l between lambda1 and lambda2 by interpolation.
	 *
	 * Returns the support of Omega that covers Reach_[t+lambda1*delta,t+lambda2*delta].
	 * To be overloaded by derived classes that can interpolate.
	 * The default implementation is to throw.
	 */
	virtual scalar Omega_support_interp_impl(const scalar& lambda1,
			const scalar& lambda2, scalar& err);

	/** Compute the support of Psi w.r.t. direction l at lambda by interpolation.
	 *
	 * Returns the support of Psi that covers Reach_[t+lambda*delta,t+lambda*delta].
	 * To be overloaded by derived classes that can interpolate.
	 * The default implementation is to throw.
	 */
	virtual scalar Psi_support_interp_impl(const scalar& lambda, scalar& err);

	/** Compute the support function of S in direction l. */
	virtual scalar
	rho(const vector_type& l, const support_function_provider& S);

	/** Compute the support vector of S in direction l.
	 *
	 * Throws if the support vector cannot be computed by S. */
	virtual vector_type
	rho_vec(const vector_type& l, const support_function_provider& S);

	const matrix_type& get_A() const {
		return handler()->my_A;
	}
	;
	const support_function_provider* get_X0() {
		return handler()->my_X0;
	}
	;
	const support_function_provider* get_U() {
		return handler()->my_U;
	}
	;
	const scalar& get_rho_X0() const {
		return my_rho_X0;
	}
	;
	const scalar& get_rho_U() const {
		return my_rho_U;
	}
	;
	const scalar& get_rho_X0_next() const {
		return my_rho_X0_next;
	}
	;

	/** Get a pointer to the handler
	 *
	 * Will be overloaded by derived classes with a pointer to the derived class.
	 * Different derived classes are not expected to work together.
	 * */
	virtual omega_model<scalar>* handler() {
		if (is_handler()) {
			return this;
		} else {
			return my_handler.get();
		}
	}

	/** Get a pointer to the handler
	 *
	 * Will be overloaded by derived classes with a pointer to the derived class.
	 * Different derived classes are not expected to work together.
	 * */
	virtual const omega_model<scalar>* handler() const {
		if (is_handler()) {
			return this;
		} else {
			return my_handler.get();
		}
	}

	/** Returns true if *this is a handler (has dynamics data) */
	bool is_handler() const;

	/** Returns true if *this has a sequence (has direction data) */
	bool has_sequence() const;

private:
	matrix_type my_A;
	scalar my_delta;
	const support_function_provider* my_X0;
	const support_function_provider* my_U;
	scalar my_rho_X0;
	scalar my_rho_U;
	scalar my_rho_X0_next;
	vector_type my_l;
	vector_type my_l_next;
	scalar my_t;

	model_ptr my_handler; // this is not an adopted pointer
	model_cache_ptr my_cache;
	scalar my_sequence_delta; // this serves to check agains the delta of the current handler to detect inconsistencies
};

template<typename scalar>
omega_model<scalar>::omega_model(const matrix_type& A, const scalar& delta,
		const support_function_provider* X0, const support_function_provider* U) :
	my_A(A), my_delta(delta), my_X0(X0), my_U(U) {
	my_handler = model_ptr(); // get_ptr() might not be available yet if it comes straight out of the factory
	// attribute it later
	my_cache = model_cache_ptr(new model_cache_type());
	//my_cache->insert(delta, get_ptr());
	my_sequence_delta = scalar(0);
}

template<typename scalar>
omega_model<scalar>::omega_model(model_ptr state_handler) :
	my_handler(state_handler), my_cache(state_handler->my_cache) {
	assert(my_handler.get() != this);
	assert(my_cache);
}

template<typename scalar>
bool omega_model<scalar>::is_handler() const {
	// handler have this pointer as null
	return !my_handler;
}

template<typename scalar>
bool omega_model<scalar>::has_sequence() const {
	// if the sequence is initialized, l is not empty
	return (my_l.size() != 0);
}

template<typename scalar>
omega_model<scalar>::~omega_model() {
}

template<typename scalar>
const typename omega_model<scalar>::matrix_type&
omega_model<scalar>::exp_AdeltaT() {
	return handler()->exp_AdeltaT_impl();
}

template<typename scalar>
scalar omega_model<scalar>::rho_X0(const vector_type& l) {
	assert(get_X0());
	return rho(l, *get_X0());
}

template<typename scalar>
scalar omega_model<scalar>::rho_U(const vector_type& l) {
	assert(get_U());
	return rho(l, *get_U());
}

template<typename scalar>
scalar omega_model<scalar>::rho(const vector_type& l,
		const support_function_provider& S) {
	scalar v;
	vector_type sv;
	bool is_empty, is_bounded;
	S.compute_support(l, v, sv, is_empty, is_bounded);
	if (!is_bounded) {
		std::stringstream s;
		logger::copyfmt_to(s);
		s << "Unbounded in direction: " << l << std::endl;
		s << "Unbounded set: " << std::endl << S;
		throw basic_exception(
				"Support function evaluation requested for an unbounded set:\n"
						+ s.str());
	}
	if (is_empty)
		throw std::runtime_error(
				"Support function evaluation requested for an empty set");
	return v;
}

template<typename scalar>
typename omega_model<scalar>::vector_type omega_model<scalar>::rho_vec(
		const vector_type& l, const support_function_provider& S) {
	scalar v;
	vector_type sv;
	bool is_empty, is_bounded;
	if (!S.computes_support_vector())
		throw std::runtime_error(
				"Support vector requested a support_function_provider that doesn't provide a support vector");
	S.compute_support(l, v, sv, is_empty, is_bounded);
	if (!is_bounded) {
		std::stringstream s;
		s << "Unbounded in direction: " << l << std::endl;
		s << "Unbounded set: " << std::endl << S;
		throw basic_exception(
				"Support function evaluation requested for an unbounded set:\n"
						+ s.str());
	}
	if (is_empty)
		throw std::runtime_error(
				"Support function evaluation requested for an empty set");
	return sv;
}

template<typename scalar>
scalar omega_model<scalar>::Omega_support(const vector_type& l, scalar& err) {
	support_sequence_ini(l);
	return support_sequence_Omega(err);
}

template<typename scalar>
void omega_model<scalar>::support_sequence_ini(const vector_type& l,
		const scalar& t0) {
	my_t = t0;
	my_l = l;
	my_l_next = exp_AdeltaT() * l;

	if (get_X0()) {
		my_rho_X0 = rho_X0(my_l);
		my_rho_X0_next = rho_X0(my_l_next);
	}
	if (get_U()) {
		my_rho_U = rho_U(my_l);
	}
	my_sequence_delta = get_delta();
	assert(my_l.size() != 0 && my_l_next.size() != 0);
}

template<typename scalar>
void omega_model<scalar>::support_sequence_next() {
	assert(my_l_next.size() != 0);
	if (has_sequence()) {
		my_t = get_t_next();
		my_l = my_l_next;
		my_l_next = exp_AdeltaT() * my_l;
		if (get_X0()) {
			my_rho_X0 = my_rho_X0_next;
			my_rho_X0_next = rho_X0(my_l_next);
		}
		if (get_U()) {
			my_rho_U = rho_U(my_l);
		}
	} else {
		throw basic_exception(
				"support_sequence_next called without initializing first");
	}
	assert(my_l.size() != 0 && my_l_next.size() != 0);
}

template<typename scalar>
scalar omega_model<scalar>::support_sequence_Omega(scalar& err) {
	if (has_sequence()) {
		return Omega_support_impl(err);
	} else {
		throw basic_exception(
				"support_sequence_Omega called without initializing first");
		return scalar();
	}
}

template<typename scalar>
scalar omega_model<scalar>::support_sequence_Psi(scalar& err) {
	if (has_sequence()) {
		return Psi_support_impl(err);
	} else {
		throw basic_exception(
				"support_sequence_Psi called without initializing first");
		return scalar();
	}
}

template<typename scalar>
bool omega_model<scalar>::support_sequence_interpolates() const {
	return false;
}

template<typename scalar>
scalar omega_model<scalar>::support_sequence_Omega_support_interp(
		const scalar& t1, const scalar& t2, scalar& err) {
	if (has_sequence()) {
		if (math::maybe(math::numeric::is_GE(t1, get_t())) && math::maybe(
				math::numeric::is_LE(t2, get_t_next())) && math::maybe(
				math::numeric::is_LE(t1, t2))) {
			scalar lambda1 = (t1 - get_t()) / get_delta();
			scalar lambda2 = (t2 - get_t()) / get_delta();

			/* Ensure lambda is between bounds */
			if (lambda1 < scalar(0))
				lambda1 = scalar(0);
			if (lambda2 > scalar(1))
				lambda2 = scalar(1);
			if (lambda1 > lambda2)
				lambda1 = lambda2;

			return Omega_support_interp_impl(lambda1, lambda2, err);
		} else {
			throw basic_exception(
					"support_sequence_Omega_support_interp called with t outside bounds");
			return scalar();
		}
	} else {
		throw basic_exception(
				"support_sequence_Omega_support_interp called without initializing first");
		return scalar();
	}
}

template<typename scalar>
scalar omega_model<scalar>::support_sequence_Psi_support_interp(
		const scalar& t, scalar& err) {
	if (has_sequence()) {
		if (math::maybe(math::numeric::is_GE(t, get_t())) && math::maybe(
				math::numeric::is_LE(t, get_t_next()))) {
			scalar lambda = (t - get_t()) / get_delta();

			/* Ensure lambda is between bounds */
			if (lambda < scalar(0))
				lambda = scalar(0);
			else if (lambda > scalar(1))
				lambda = scalar(1);

			return Psi_support_interp_impl(lambda, err);
		} else {
			throw basic_exception(
					"support_sequence_Psi_support_interp called with t outside bounds");
			return scalar();
		}
	} else {
		throw basic_exception(
				"support_sequence_Psi_support_interp called without initializing first");
		return scalar();
	}
}

template<typename scalar>
scalar omega_model<scalar>::Omega_support_interp_impl(const scalar& lambda1,
		const scalar& lambda2, scalar& err) {
	throw basic_exception(
			"support_sequence_Omega_support_interp called but model doesn't provide it");
}

template<typename scalar>
scalar omega_model<scalar>::Psi_support_interp_impl(const scalar& lambda,
		scalar& err) {
	throw basic_exception(
			"support_sequence_Psi_support_interp called but model doesn't provide it");
}
template<typename scalar>
scalar omega_model<scalar>::Psi_support(const vector_type& l, scalar& err) {
	support_sequence_ini(l);
	return support_sequence_Psi(err);
}

template<typename scalar>
void omega_model<scalar>::update(const scalar& delta) {
	// Step 1: define delta for the dynamics memory

	// if *this is a handler, update
	if (is_handler()) {
		if (!math::numeric::is_MEQ(my_delta, delta)) {
			my_delta = delta;
			update_impl();
		}
	} else {
		// if *this is not a handler, go find a handler
		if (!math::numeric::is_MEQ(get_delta(), delta)) {
			typename model_cache_type::const_iterator it =
					my_cache->find(delta);
			if (it != my_cache->end()) {
				my_handler = it->second;
			} else {
				// create a copy of the current handler, update it, and add it to the cache
				my_handler = model_ptr(my_handler->clone());
				my_handler->update(delta);
				my_cache->insert(delta, my_handler);
			}
		}
	}

	// Step 1: define delta for the sequence memory

	// update the next state in the sequence
	if (has_sequence()) {
		if (!math::numeric::is_MEQ(my_sequence_delta, delta)) {
			my_l_next = exp_AdeltaT() * my_l;
			if (get_X0()) {
				my_rho_X0_next = rho_X0(my_l_next);
			}
			my_sequence_delta = delta;
		}

		// consistency check
		if (!math::numeric::is_MEQ(my_sequence_delta, get_delta())) {
			throw basic_exception(
					"handler delta and sequence delta are different");
		}
	}


}

template<typename scalar>
const typename omega_model<scalar>::vector_type& omega_model<scalar>::get_l() const {
	return my_l;
}

template<typename scalar>
const typename omega_model<scalar>::vector_type& omega_model<scalar>::get_l_next() const {
	return my_l_next;
}

template<typename scalar>
const scalar& omega_model<scalar>::get_delta() const {
	return handler()->my_delta;
}

template<typename scalar>
const scalar& omega_model<scalar>::get_t() const {
	return my_t;
}

template<typename scalar>
scalar omega_model<scalar>::get_t_next() const {
	return my_t + my_sequence_delta;
}

template<typename scalar>
typename omega_model<scalar>::model_ptr omega_model<scalar>::get_ptr() {
	model_ptr
			p =
					boost::enable_shared_from_this<omega_model<scalar> >::shared_from_this();
	return p;
}

}
}

#endif /* OMEGA_MODEL_H_ */
