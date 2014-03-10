/*
 * sf_evaluator.hpp
 *
 *  Created on: Jan 20, 2011
 *      Author: frehse
 */

#ifndef SF_EVALUATOR_HPP_
#define SF_EVALUATOR_HPP_

#include <sstream>
#include "utility/basic_exception.h"
#include "utility/basic_warning.h"
#include "utility/logger_stopwatch.h"

#include "math/numeric/container_comp.h"
#include "core/continuous/support_function/sf_base/sf_unary.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/continuous/support_function/sf_derived/sf_ode_deterministic.hpp"
#include "core/continuous/support_function_provider_utility.h"

namespace continuous {
namespace support_function {

template<class T>
sf_evaluator<T>::sf_evaluator(math::matrix<T> nA, math::vector<T> nb,
		index_to_variable_id_map_ptr niimap,
		support_function_provider::const_ptr nX0,
		support_function_provider::const_ptr nU) :
	A(nA), b(nb), X0(nX0), U_orig(nU), dom(niimap) {
	assert(A.size1() == A.size2());
	assert(A.size2() == b.size());
	assert(A.size1() == dom.size());
	assert(!math::definitely(X0->is_empty()));

	Am = typename omega_model<T>::matrix_type(dom, dom, A);

	use_simple_err = true;

	delta = T(0);

	U = U_orig;

	if (!U || math::definitely(U->is_empty())) {
		// There are no inputs, so U = {b} and centered U = {0}
		math::vector<T> zeros(dom.size(), T(0));
		U = support_function_provider::const_ptr(new finite_hyperbox<T> (b, zeros, dom));
/*		U_centered = support_function_provider::const_ptr(new hyperbox<T> (
				zeros, zeros, iimap));
*/
		U_centered = support_function_provider::const_ptr();
		u0 = b;
	} else {
		math::vdom_vector<T> b_named(dom, b);

		// add b to U
		if (!math::numeric::is_MEQ(b, T(0))) {
			math::affine_map<T> translate_map(dom, b_named);
			U = support_function_provider::const_ptr(
					new support_function::sf_unary<T>(U, translate_map));
		}

		finite_hyperbox<T> U_box;
		typename finite_hyperbox<T>::vdom_vector_type U_center_tmp;
		try {
			U_box = finite_bounding_box<T> (*U);
			U_center_tmp = U_box.get_c_dom();
		} catch (std::exception& e) {
			std::stringstream ss;
			logger::copyfmt_to(ss);
//			ss << *U;
			throw basic_exception("Can't compute center of the set of values u in dynamics x'==Ax+u. Did you forget to bound the input variables in the invariant?\n", e);
		}
		U_center_tmp.remap(dom);

		// obtain centered U
		if (math::numeric::is_MEQ(U_center_tmp.get_vector(), b)) {
			// U was already centered if the center is b
			U_centered = U_orig;
		} else if (!math::numeric::is_MEQ(U_center_tmp.get_vector(), T(0))) {
			// If nonzero center then shift.
			// In order not to cascade two sf_unary shifts,
			// do the shift on the original.
			math::affine_map<T> translate_map(dom, b_named - U_center_tmp);

			U_centered = support_function_provider::const_ptr(
					new support_function::sf_unary<T>(U_orig, translate_map));

//			std::cout << "translate map: "<< translate_map << std::endl;
//			finite_hyperbox<T> dummy = finite_bounding_box<T> (*U_orig);
//			std::cout << "U_orig bounding box: " << dummy << std::endl;
//			dummy = finite_bounding_box<T> (*U_centered);
//						std::cout << "U_orig centered bounding box: " << dummy << std::endl;
		} else {
			U_centered = U;
		}

		////		std::cout << "U_center: "<< U_center_tmp << std::endl;
		//		finite_hyperbox<T> dummy = finite_bounding_box<T> (*U_centered);
		//		std::cout << "U_centered bounding box: " << dummy << std::endl;
		//		dummy = finite_bounding_box<T> (*U);
		//		std::cout << "U bounding box: " << dummy << std::endl;

		T mu_centered;
		get_max_infinity_norm(*U_centered, mu_centered);
		if (math::numeric::is_MEQ(mu_centered, T(0))) {
			U_centered = support_function_provider::const_ptr();
		}
		u0 = U_center_tmp.get_vector();
	}
}

template<class T>
T sf_evaluator<T>::get_timestep(unsigned int i, T X_bloat_value,
		T V_bloat_value, double step_tolerance) {
	using namespace math;
	using namespace math::numeric;

	const double fac = 2.0; //1.5; // 1.05; //

	T d = delta;
	T tol(step_tolerance);
	if (step_tolerance > 0.0) {
		T error_val = X_bloat_value + V_bloat_value;
		if (definitely(is_GT(error_val, tol)))
			d = delta / T(fac);
		else if (definitely(is_LT(
				T(fac) * T(fac) * T(fac) * T(fac) * error_val, tol))) {
			//d = T(fac) * T(fac) * delta;
			d = T(fac) * delta;
		}
	}

	//	static double fac_exp = 0;
	//
	//	T d = delta;
	//	T tol(step_tolerance);
	//	if (step_tolerance > 0.0) {
	//		T error_val = std::max(X_bloat_value, V_bloat_value);
	//		double fac = std::pow((1.0+std::abs(fac_exp)/10.0),fac_exp);
	//		double fac_plus = std::pow((1.0+std::abs(fac_exp+1)/10.0),fac_exp+1);
	//		if (definitely(is_GT(error_val, tol)))
	//			fac_exp=fac_exp-1.0;
	//		else if (definitely(is_LT(T(fac_plus/fac)*error_val, tol))) {
	//			fac_exp=fac_exp+1.0;
	//		}
	//		fac = std::pow((1.0+std::abs(fac_exp)/10.0),fac_exp);
	//		std::cout << fac << " exp " << fac_exp << std::endl;
	//		d = orig_delta * T(fac);
	//	}
	return d;
}

template<class T>
math::matrix<T> sf_evaluator<T>::compute_matrix_tol(const std::list<
		math::vector<T> >& directions, const std::vector<
		scalar_with_infinity<T> > inh_coeffs, const T& ndelta,
		const T& time_horizon, std::vector<T>& delta_vec, double step_tolerance) {
	assert(directions.empty() || A.size1() == directions.begin()->size());
	assert(inh_coeffs.size() == directions.size());

	LOGGERSWOC(DEBUG4,"sf_evaluator<T>::compute_matrix_tol","Computing support function sequence with tolerance");

	set_delta(ndelta);
	delta_vec = std::vector<T>();

	unsigned int R = directions.size();

	// define tolerances for omega and psi,
	// i.e., distribute error on initial states and inputs
	T omega_tolerance;
	T psi_tolerance;
	if (U) {
		omega_tolerance = T(0.5)*step_tolerance;
		psi_tolerance = T(0.5)*step_tolerance;
	} else {
		// no inputs, everyhting is for omega
		omega_tolerance = step_tolerance;
		psi_tolerance = T(0);
	}

	// loop over MAX_ITERS and break when an \omega_i set goes completely out in any of the invariant direction for the first time
	LOGGER(DEBUG3,"sf_evaluator<T>::compute_matrix","loop over iterations");
	//		std::cout << "Inside loop over MAX_ITERS" << std::endl;
	std::list<std::vector<T> > my_sfm_cols;
	std::vector<T> si(R, T(0));
	std::vector<T> ci(R, T(0));
	std::vector<T> rho_Omegai(R, T(0));
	//std::vector<T> delta_sugg(R, delta);
	std::vector<T> Omega_bloat_value(R, T(0));
	std::vector<T> Psi_bloat_value(R, T(0));
	// accumulation error for each direction
	math::vector<T> err_Psi(R,T(0));

	// Approximation models for each direction, R being the number of directions
	std::vector<omega_model_ptr> Omega_models(R);
	std::vector<omega_model_ptr> Psi_models(R);

	//ri's initially contains the original directions.

	// Create the state model from which all others are cloned
	{
		/** This is wrapped in a while loop because create_models could crash */
		bool success = false;
		while (!success) {
			try {
				create_models(delta);
				success = true;
			} catch (math::invalid_number_exception& e) {
				set_delta(delta / 16);
			}
		}
	}

	{
		unsigned int j = 0;
		for (typename std::list<math::vector<T> >::const_iterator it =
				directions.begin(); it != directions.end(); ++it, j++) {
			Omega_models[j] = omega_model_ptr(Omega_model->client());
			Psi_models[j] = omega_model_ptr(Psi_model->client());
			Omega_models[j]->support_sequence_ini(
					math::vdom_vector<T>(dom, *it));
			Psi_models[j]->support_sequence_ini(math::vdom_vector<T>(dom, *it));
		}
	}

	T global_bloat_max(0); // globally largest error
	unsigned int global_bloat_max_dir(0); // direction with largest error

	T time = T(0);
	T delta_min;

	std::vector<T> old_si(si);
	std::vector<T> old_ci(ci);

	bool inv_violation = false;
	unsigned int i = 0;
	while (!inv_violation && time < time_horizon) {

		T bloat_max_val = T(0);
		unsigned int bloat_max_dir(0); // direction with largest error

		delta_min = T(999) * delta;
		for (unsigned int j = 0; !inv_violation && j < directions.size(); j++) {
			T rho_Omega0 = Omega_models[j]->support_sequence_Omega(
					Omega_bloat_value[j]); // \rho_{Omega_{0}}{r}
			rho_Omegai[j] = rho_Omega0 + si[j] + ci[j];

//			std::cout << "in domain " << dom << " current l(" << j << ")"
//					<< Omega_models[j]->get_l() << ", rho_Omega0="
//					<< rho_Omega0 << ", rho_Psi=" << si[j] << ", rho_Ucenter="
//					<< ci[j] << std::endl;

			if (!inh_coeffs[j].is_pos_infinity() && math::definitely(
					math::numeric::is_GT(-rho_Omegai[j],
							inh_coeffs[j].get_val()))) {
				inv_violation = true;
				IFLOGGER(DEBUG4) {
					std::stringstream ss;
					logger::copyfmt_to(ss);
					// get the current direction
					typename std::list<math::vector<T> >::const_iterator lit;
					unsigned int lj = 0;
					for (lit = directions.begin(); lit != directions.end()
							&& lj < j; ++lit) {
						++lj;
					}
					ss << math::lin_expression<T>(dom, -*lit, T(0)) << " > " << inh_coeffs[j].get_val()
						 << " with min surpassing " << -rho_Omegai[j];
					LOGGER(DEBUG4,"compute_matrix","at step "+to_string(i)+", t="+to_string(time)+", violated inv because "+ss.str());
				}
			}
			T rho_Psi0 = Psi_models[j]->support_sequence_Psi(Psi_bloat_value[j]); // \rho_{V}{r}

//			std::cout << "in domain " << dom << " current l("<<j<<")" << Psi_models[j]->get_l()
//					<< " with center shift " << math::vdom_vector<T>(dom,
//					U_center_shifted) << ": adding " << scalar_product(
//					Psi_models[j]->get_l(),
//					math::vdom_vector<T>(dom, U_center_shifted)) << std::endl;

			ci[j] = ci[j] + scalar_product(Psi_models[j]->get_l(),
					math::vdom_vector<T>(dom, U_center_shifted));
			si[j] = si[j] + rho_Psi0;

			// update time step
			T err_margin = omega_tolerance + psi_tolerance * (time + delta)
					/ time_horizon;
			// need to take into account the psi error we'll incur during the current step
			// (even if it doesn't figure until the next step)
			T delta_sugg = get_timestep(i, Omega_bloat_value[j], err_Psi[j] + Psi_bloat_value[j],
					err_margin);
			if (delta_sugg < delta_min)
				delta_min = delta_sugg;

			T error_val = Omega_bloat_value[j] + err_Psi[j];
			if (error_val > bloat_max_val) {
				bloat_max_val = error_val;
				bloat_max_dir = j;
			}
			bloat_max_val = std::max(bloat_max_val, error_val);
			//std::cout << "nX:" << Omega_bloat_value[j] << "nV:" << V_bloat_value << " e:" << error_val << " max: " << bloat_max_val << std::endl;
		}
		if (!inv_violation) {
			if (math::maybe(math::numeric::is_GE(delta_min, delta))) {
				// keep track of worst error
				if (bloat_max_val > global_bloat_max) {
					global_bloat_max = bloat_max_val;
					global_bloat_max_dir = bloat_max_dir;
				}

				// accept the new values
				time += delta;
				delta_vec.push_back(delta);
				my_sfm_cols.push_back(rho_Omegai);

				for (unsigned int jj = 0; jj < directions.size(); jj++) {
					err_Psi[jj] += Psi_bloat_value[jj];
					Omega_models[jj]->support_sequence_next();
					Psi_models[jj]->support_sequence_next();
				}

				old_si=si;
				old_ci=ci;

				++i;

				// check if we can increase the time values
				if (math::definitely(math::numeric::is_GT(delta_min, delta))) {
					LOGGER(DEBUG3,"compute_matrix","at step "+to_string(i)+", t="+to_string(time)+", err.<="+to_string(bloat_max_val)+", +time step:="+to_string(delta_min));
					set_delta(delta_min); // for U_center_shifted
					for (unsigned int jj = 0; jj < directions.size(); jj++) {
						set_delta(*Omega_models[jj], delta);
						set_delta(*Psi_models[jj], delta);
					}
				}
			} else {
				LOGGER(DEBUG3,"compute_matrix","at step "+to_string(i)+", t="+to_string(time)+", err.<="+to_string(bloat_max_val)+", -time step:="+to_string(delta_min));
				// redo this turn
				si = old_si;
				ci = old_ci;
				set_delta(delta_min); // for U_center_shifted
				for (unsigned int jj = 0; jj < directions.size(); jj++) {
					set_delta(*Omega_models[jj], delta);
					set_delta(*Psi_models[jj], delta);
				}
				// we diminished the time step so we don't advance to the next
				// item in the support_sequence
			}
		}
	}
	// show the error and direction with this error
	IFLOGGER(DEBUG2) {
		std::stringstream ss;
		logger::copyfmt_to(ss);
		// get the currenct direction
		typename std::list<math::vector<T> >::const_iterator lit;
		unsigned int lj;
		for (lit = directions.begin(); lit != directions.end()
				&& lj < global_bloat_max_dir; ++lit, lj++) {
		}
		ss << math::lin_expression<T>(dom, *lit, T(0));
		LOGGER(DEBUG2,"compute_matrix","max error "+to_string(global_bloat_max)+ " in direction "+ss.str());
	}


	unsigned int cols = my_sfm_cols.size();
	math::matrix<T> Y(R, cols);

	// Constructing the SFM
	unsigned int col = 0;
	for (typename std::list<std::vector<T> >::const_iterator it =
			my_sfm_cols.begin(); it != my_sfm_cols.end(); ++it, col++) {
		for (unsigned int i = 0; i < R; i++)
			Y(i, col) = (*it)[i];
	}

	// We're supposed to compute time elapse until the invariant is violated.
	// If it's not violated, this means we hit a timeout.
	// Let's warn the user that the result is incomplete
	if (!inv_violation) {
		basic_warning(
				"get_support_function_evaluations",
				"Reached time horizon without exhausting all states, result is incomplete.",
				basic_warning::INCOMPLETE_OUTPUT);
	}
	return Y;
}

template<class T>
math::matrix<T> sf_evaluator<T>::compute_matrix(
		const std::list<math::vector<T> >& directions, const T& ndelta,
		unsigned int N, std::vector<T> delta_vec) {
	assert(directions.empty() || A.size1() == directions.begin()->size());
	// @todo this isn't yet so since we don't yet use time stamps but deltas assert(delta_vec.size() == 0 || delta_vec.size() >= N + 1);

	LOGGERSWOC(DEBUG4,"sf_evaluator<T>::compute_matrix_tol","Computing support function sequence for already existing time points");

	bool use_delta_vec = delta_vec.size() > 0;
	set_delta(ndelta);
	create_models(delta);

	unsigned int R = directions.size();
	math::vector<T> r; // direction
	T s; // input sf value

	unsigned int row = 0;
	T rho_Omega0 = T(0);
	T rho_Psi0 = T(0);
	T X_bloat_value = T(-1.0);
	T V_bloat_value = T(-1.0);
	T bloat_max_val = T(0);

	// for computing the sequence OSum_i \Phi^i U_center_shifted
	T c(0);

	logger::logger_id
			logid =
					LOGGER(DEBUG6,"sf_evaluator<T>::compute_matrix","fixed horizon support matrix computation ");

	//		std::cout << "Inside loop over directions\n" << std::endl;
	// main loop ranging over directions
	math::matrix<T> Y(R, N); // Y initialised to a 2-dim matrix of size R X N, R being the number of directions.

	for (typename std::list<math::vector<T> >::const_iterator it =
			directions.begin(); it != directions.end(); it++, row++) {
		r = *it;

		Omega_model->support_sequence_ini(math::vdom_vector<T>(dom, r));
		Psi_model->support_sequence_ini(math::vdom_vector<T>(dom, r));
		if (use_delta_vec) {
			set_delta_unchecked(delta_vec[0]); // for U_center_shifted
//			set_delta(*Omega_model, delta_vec[0]); // compute all matrices
//			set_delta(*Psi_model, delta_vec[0]); // compute all matrices
		}
		//r = math::vector<T>(*it); // r initialised to l
		//
		s = T(0); // s initialised to 0
		c = T(0); // c initialised to 0

		for (unsigned int i = 0; i < N; i++) {
			rho_Omega0 = Omega_model->support_sequence_Omega(X_bloat_value);

			Y(row, i) = rho_Omega0 + s + c;
			rho_Psi0 = Psi_model->support_sequence_Psi(V_bloat_value); // \rho_{V}{r}
			c = c + scalar_product(Psi_model->get_l(), math::vdom_vector<T>(
					dom, U_center_shifted));
			s = s + rho_Psi0;


			bloat_max_val = std::max(bloat_max_val, X_bloat_value
					+ V_bloat_value);
			if (i + 1 < N) {
				Omega_model->support_sequence_next();
				Psi_model->support_sequence_next();
				if (use_delta_vec && !math::numeric::is_MEQ(delta_vec[i],
						delta_vec[i + 1])) {
//std::cout << "dir: "<<*it << " iter "<<i << " delta " << delta_vec[i + 1] << Omega_model->get_delta() << ", rho_omega0: " << rho_Omega0 << ", rho_V: " << s+c << std::endl;
					set_delta_unchecked(delta_vec[i + 1]); // for U_center_shifted
//					set_delta(*Omega_model, delta_vec[i + 1]); // compute all matrices
//					set_delta(*Psi_model, delta_vec[i + 1]); // compute all matrices
				}
			}
		}
	}

	LOGGER_ATTACH(DEBUG6,"compute_matrix","finished with err.<="+to_string(bloat_max_val),logid);

	return Y;
}

template<class T>
void sf_evaluator<T>::set_delta(omega_model<T>& omega_m, T d) {
	omega_m.update(d);
}

template<class T>
void sf_evaluator<T>::set_delta(T d) {
	bool success = false;
	while (!success) {
		try {
			set_delta_unchecked(d);
			success = true;
		} catch (math::invalid_number_exception& e) {
			d = d / T(16.0);
		}
	}
}

template<class T>
void sf_evaluator<T>::set_delta_unchecked(T d) {
	delta = d;
	/* New computation of Phi at the same time as
	 * special matrices
	 * Phi=e^(A*delta)
	 * Phi1=A^{-1}(e^(A*delta)-I)
	 * Phi2=A^{-2}(e^(A*delta)-I-delta*A)
	 */
	typename matrix_cache_type::const_iterator it =
			my_cache.find(delta);
	if (it != my_cache.end()) {
		Phi1 = it->second;
	} else {
		math::get_special_matrices(A, delta, Phi, Phi1, Phi2);
		my_cache.insert(delta, Phi1);
	}
	// @note The Phi1 formula for shifting the center is valid even for zero A matrix
	//	if (math::numeric::is_MEQ(A,T(0))) {
//		U_center_shifted = delta * u0;
//	} else
		U_center_shifted = Phi1 * u0;


	// update models
	if (Omega_model)
		Omega_model->update(delta);
	if (Psi_model)
			Psi_model->update(delta);
}

template<class T>
void sf_evaluator<T>::create_models(const T& delta) {
	if (!Omega_model) {
		Omega_model = omega_model_ptr(
				omega_model_factory::create(Am, delta, X0.get(), U.get()));
	} else {
		Omega_model->update(delta);
	}
	if (!Psi_model) {
		Psi_model = omega_model_ptr(
				omega_model_factory::create(Am, delta, 0, U_centered.get()));
	} else {
		Psi_model->update(delta);
	}
}

}
}


#endif /* SF_EVALUATOR_HPP_ */
