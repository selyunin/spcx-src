/*
 * sf_evaluator.h
 *
 *  Created on: Jan 20, 2011
 *      Author: frehse
 */

#ifndef SF_EVALUATOR_H_
#define SF_EVALUATOR_H_

#include "math/matrix.h"
#include "math/matrix_operators.h"
#include "math/vector.h"
#include "math/scalar_types/scalar_with_infinity.h"
#include "core/continuous/polyhedra/polyhedron.h"
#include "core/continuous/polyhedra/polyhedron_utility.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/ppl_polyhedron/convert_to_ppl.h"
#include "core/discrete/discrete_set_stl_set.h"
#include "core/continuous/polyhedra/hyperbox/bounding_box.h"
#include "core/continuous/support_function/template_directions/choose_directions.h"

#include "omega_model_factory.h"

namespace continuous {
namespace support_function {

template<class T>
class sf_evaluator {

public:
	/** Construct an evaluator with A,b,X0,U,iimap. */
	sf_evaluator(math::matrix<T> nA, math::vector<T> nb,
			index_to_variable_id_map_ptr niimap,
			support_function_provider::const_ptr nX0,
			support_function_provider::const_ptr nU);

	/**
	 * \brief Evaluate the support function in the directions given by \p directions,
	 * for N time steps separated by delta time units each. The initial set is X0,
	 * the dynamics are given by the matrix A. The set U defines the bounds on the inputs
	 * (in continuous time).
	 *
	 * @param inh_coeffs A vector to denote which directions are invariant directions. The
	 * entries i which are other than pos_infty implies that the corresponding directions i
	 * in the directions list are invariant directions.
	 * For every invariant direction l with corresponding inh_coeff c, the states must satisfy
	 * l^Tx <= c.
	 *
	 * delta_vec[i] is the time duration for each set Omega_i, i from 0,...,N-1.
	 */
	math::matrix<T> compute_matrix_tol(
			const std::list<math::vector<T> >& directions, const std::vector<
					scalar_with_infinity<T> > inh_coeffs, const T& delta,
			const T& time_horizon, std::vector<T>& delta_vec,
			double step_tolerance = -1.0);

	/**
	 * \brief Evaluate the support function in the directions given by \p directions,
	 * for N time steps separated by delta time units each. The initial set is X0,
	 * the dynamics are given by the matrix A. The set U defines the bounds on the inputs
	 * (in continuous time).
	 *
	 * @param delta_vec If provided, must be of at least size N-1 and specifies the delta for
	 * each time step from 0,...,n-1.
	 */
	math::matrix<T> compute_matrix(
			const std::list<math::vector<T> >& directions, const T& delta,
			unsigned int N, std::vector<T> delta_vec = std::vector<T>());

private:
	bool use_simple_err;

	typedef math::matrix<T> matrix_type;
	typedef math::vector<T> vector_type;
	typedef boost::shared_ptr<omega_model<T> > omega_model_ptr;

	// Parameters of the algorithm
	matrix_type A;
	vector_type b;
	support_function_provider::const_ptr X0;
	support_function_provider::const_ptr U;
	support_function_provider::const_ptr U_orig;
	support_function_provider::const_ptr U_centered;
	vector_type u0; // the center of U
	positional_vdomain dom;

	// Internal parameters that depend on delta
	T delta;
	matrix_type Phi; //=e^(A*delta)
	matrix_type Phi_T; //=Phi^T
	matrix_type Phi1; //=A^{-1}(e^(A*delta)-I)
	matrix_type Phi2; //=A^{-2}(e^(A*delta)-I-delta*A)
	vector_type U_center_shifted;

	/** Set delta and compute parameters depending on delta with numerical checks.
	 *
	 * If numerical problems arise, the process is automatically
	 * repeated with a smaller delta.
	 * @attention  It holds that afterwards delta<=d, but it
	 * may be that delta is much smaller than d for numerical reasons. */
	void set_delta(T d);

	/** Set delta and compute parameters depending on delta without numerical checks.
	 *
	 * @attention Only to be used by set_delta. Do not use directly. */
	void set_delta_unchecked(T d);

	/** Set delta in the omega_model. */
	void set_delta(omega_model<T>& omega_m, T d);

	/** Adapt the timestep at iteration i according to the bloating factor
	 * of Omega_0 and the support value of V.
	 */
	T get_timestep(unsigned int i, T bloat_value, T V_support_value,
			double step_tolerance);

	/** Create omega and psi models */
	void create_models(const T& delta);

	// For caching Phi1
	typedef math::unique_scalar_to_value_store<T, matrix_type>
			matrix_cache_type;
	matrix_cache_type my_cache;

	// For time step estimate
	T step_base;
	T last_error;
	T last_delta;
	double fac;

	// For omega and psi models */
	typename omega_model<T>::matrix_type Am;
	omega_model_ptr Omega_model;
	omega_model_ptr Psi_model;
};

}
}

#include "sf_evaluator.hpp"

#endif /* SF_EVALUATOR_H_ */
