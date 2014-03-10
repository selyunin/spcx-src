/*
 * ode_solver_cvode.h
 *
 *  Created on: Oct 7, 2010
 *      Author: frehse
 */

#ifndef ODE_SOLVER_CVODE_H_
#define ODE_SOLVER_CVODE_H_

#include <boost/shared_ptr.hpp>
#include "math/vdom/trajectory.h"
#include "math/vdom/state_functor_utility.h"
#include "math/vdom/state_vector_functor_utility.h"
#include "math/vdom/state_matrix_functor_utility.h"
#include "ode_solver.h"

#include <sundials/sundials_config.h> /* definition of type realtype */

namespace math {
namespace ode {

/** Mimic the typedef from cvode without having to include the sundial_types.h
 * header (which would pollute our namespace).
 */

namespace sundial {
#if defined(SUNDIALS_SINGLE_PRECISION)

typedef float realtype;
# define RCONST(x) x##F
# define BIG_REAL FLT_MAX
# define SMALL_REAL FLT_MIN
# define UNIT_ROUNDOFF FLT_EPSILON

#elif defined(SUNDIALS_DOUBLE_PRECISION)

typedef double realtype;
# define RCONST(x) x
# define BIG_REAL DBL_MAX
# define SMALL_REAL DBL_MIN
# define UNIT_ROUNDOFF DBL_EPSILON

#elif defined(SUNDIALS_EXTENDED_PRECISION)

typedef long double realtype;
# define RCONST(x) x##L
# define BIG_REAL LDBL_MAX
# define SMALL_REAL LDBL_MIN
# define UNIT_ROUNDOFF LDBL_EPSILON

#endif
}

/** An ode solver interface to cvode.
 *
 * cvode is available at http://computation.llnl.gov/casc/sundials/
 */
template<typename scalar_type> class ode_solver_cvode: public ode_solver<
		scalar_type> {
public:
	typedef boost::shared_ptr<ode_solver<scalar_type> > ptr;
	typedef boost::shared_ptr<const ode_solver<scalar_type> > const_ptr;

	typedef typename ode_solver<scalar_type>::state state;
	typedef typename math::state_functor<scalar_type> state_functor;
	typedef typename math::state_vector_functor<scalar_type> state_vector_functor;
	typedef typename math::state_matrix_functor<scalar_type> state_matrix_functor;
	typedef typename ode_solver<scalar_type>::ode_parameters ode_parameters;
	typedef typename ode_solver<scalar_type>::ivp_result ivp_result;
	typedef typename ode_solver<scalar_type>::rootf_result rootf_result;

	typedef sundial::realtype implementor_scalar;
	typedef ode_solver_cvode<implementor_scalar> implementor_type;

	ode_solver_cvode();
	~ode_solver_cvode();

	/** Solve the currently defined initial value problem up to time tfinal.
	 * May return just the final state, even when tfinal>max_timestep.
	 */
	virtual ivp_result step_ivp(scalar_type tfinal);

	/** Solve the currently defined root finding problem up to time tfinal,
	 * stopping as soon as a root is found.
	 *
	 * May return just the final state, even when tfinal>max_timestep.
	 */
	virtual rootf_result step_rootf(scalar_type tfinal);

	/**
	 * changes the root function of the current ivp to g
	 * WARNING :this function should not be called prior to a call
	 * to init_rootf.
	 * @param g, the new rootfinding function
	 */
	void change_rootf(const state_vector_functor& g);

	/**
	 * Used to specify a functor used to produce the Jacobian
	 * matrix for the solver (instead of an internal approximation)
	 */
	void specify_jacobian(const state_matrix_functor & j);

	virtual scalar_type get_solver_current_time() const;

	/** Get a string name for the solver. */
	virtual std::string get_solver_name() const {
		return "CVODE";
	}
	;

protected:
	virtual void init_ivp_impl();
	virtual void init_rootf_impl();

private:
	implementor_type* my_implementor;
	typename ode_solver<implementor_scalar>::state_functor* my_conv_f;
	typename ode_solver<implementor_scalar>::state_vector_functor* my_conv_g;
	typename ode_solver<implementor_scalar>::state_matrix_functor* my_conv_j;
};

template<>
class ode_solver_cvode<sundial::realtype> : public ode_solver<sundial::realtype> {
public:
	typedef boost::shared_ptr<ode_solver_cvode<sundial::realtype> > ptr;
	typedef boost::shared_ptr<const ode_solver_cvode<sundial::realtype> >
			const_ptr;

	typedef ode_solver<sundial::realtype>::state state;
	typedef ode_solver<sundial::realtype>::state_functor state_functor;
	typedef ode_solver<sundial::realtype>::state_vector_functor state_vector_functor;
	typedef ode_solver<sundial::realtype>::vector vector;

	struct impl_data;

	ode_solver_cvode();
	~ode_solver_cvode();

	/** Solve the currently defined initial value problem up to time tfinal.
	 */
	virtual ivp_result step_ivp(sundial::realtype tfinal);

	/** Solve the currently defined root finding problem up to time tfinal.
	 */
	virtual rootf_result step_rootf(sundial::realtype tfinal);

	/**
	 * changes the root function of the current ivp to g
	 * WARNING :this function should not be called prior to a call
	 * to init_rootf
	 */
	void change_rootf(const state_vector_functor& g);

	/**
	 * Used to specify a functor used to produce the Jacobian
	 * matrix for the solver (instead of an internal approximation)
	 */
	virtual void specify_jacobian(const state_matrix_functor & j);

	virtual sundial::realtype get_solver_current_time() const;
	/** Get a string name for the solver. */
	virtual std::string get_solver_name() const;

protected:
	virtual void init_ivp_impl();
	virtual void init_rootf_impl();

private:
	impl_data* my_data;
};

/** Implementations of the ode_solver_cvode.
 * They are attached at the end to avoid instantiation of the implementor class
 * before the specializations are defined.
 */
template<typename scalar_type>
ode_solver_cvode<scalar_type>::ode_solver_cvode() :
	my_implementor(new implementor_type()), my_conv_f(0), my_conv_g(0), my_conv_j(0) {
}

template<typename scalar_type>
ode_solver_cvode<scalar_type>::~ode_solver_cvode() {
	delete my_implementor;
	delete my_conv_f;
	delete my_conv_g;
	delete my_conv_j;
}

template<typename scalar_type>
void ode_solver_cvode<scalar_type>::init_ivp_impl() {
	typename implementor_type::state x0impl = this->my_x0.template convert_to<
			implementor_scalar> ();
	implementor_scalar t0impl = convert_element<implementor_scalar> (
			this->my_t0);
	delete my_conv_f;
	delete my_conv_g;
	my_conv_f = new convert_state_functor<scalar_type,implementor_scalar>(
			*this->my_f);
	typename implementor_type::ode_parameters par(this->my_odepar.template convert_to<implementor_scalar>());
	my_implementor->init_ivp(x0impl, t0impl, *this->my_conv_f, par);
}

template<typename scalar_type>
typename ode_solver_cvode<scalar_type>::ivp_result ode_solver_cvode<scalar_type>::step_ivp(
		scalar_type tfinal) {
	implementor_scalar tfinalimpl =
			convert_element<implementor_scalar> (tfinal);
	typename implementor_type::ivp_result resimpl = my_implementor->solve_ivp(
			tfinalimpl);
	ivp_result res = resimpl.template convert_to<scalar_type> ();
	return res;
}

template<typename scalar_type>
void ode_solver_cvode<scalar_type>::init_rootf_impl() {
	typename implementor_type::state x0impl = this->my_x0.template convert_to<
			implementor_scalar> ();
	implementor_scalar t0impl = convert_element<implementor_scalar> (
			this->my_t0);
	delete my_conv_f;
	delete my_conv_g;
	my_conv_f = new convert_state_functor<scalar_type,implementor_scalar>(
			*this->my_f);
	my_conv_g = new convert_state_vector_functor<scalar_type,implementor_scalar>(
			*this->my_g);
	typename implementor_type::ode_parameters par(this->my_odepar.template convert_to<implementor_scalar>());
	my_implementor->init_rootf(x0impl, t0impl, *this->my_conv_f,
			*this->my_conv_g, par);
}

template<typename scalar_type> typename ode_solver_cvode<scalar_type>::rootf_result ode_solver_cvode<
		scalar_type>::step_rootf(scalar_type tfinal) {
	implementor_scalar tfinalimpl =
			convert_element<implementor_scalar> (tfinal);
	typename implementor_type::rootf_result resimpl =
			my_implementor->step_rootf(tfinalimpl);
	rootf_result res = resimpl.template convert_to<scalar_type> ();
	return res;
}

template<typename scalar_type>
void ode_solver_cvode<scalar_type>::change_rootf(const state_vector_functor & g) {

	delete my_conv_g;
	my_conv_g = new convert_state_vector_functor<scalar_type,implementor_scalar>(
			*this->my_g);
	my_implementor->change_rootf(*this->my_conv_g);

}
template<typename scalar_type>
void ode_solver_cvode<scalar_type>::specify_jacobian(const state_matrix_functor & j){
	delete my_conv_j;
	my_conv_j = new convert_state_matrix_functor<scalar_type,implementor_scalar>(
			*this->my_j);
	my_implementor->specify_jacobian(*this->my_conv_j);

}

template <typename scalar_type>
scalar_type ode_solver_cvode<scalar_type>::get_solver_current_time() const
{
		return convert_element<scalar_type> (my_implementor->get_solver_current_time());
		//return typename scalar_type::convert_to<typename implementor_type>(my_implementor->get_solver_current_time());
}


}
}

#endif /* ODE_SOLVER_CVODE_H_ */
