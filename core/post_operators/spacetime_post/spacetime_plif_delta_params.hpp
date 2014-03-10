/*
 * spacetime_plif_delta_params.hpp
 *
 *  Created on: Oct 20, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_PLIF_DELTA_PARAMS_HPP_
#define SPACETIME_PLIF_DELTA_PARAMS_HPP_

#include "spacetime_plif.h"

#include "core/continuous/support_function/sf_base/sf_unary.h"

namespace spacetime {

inline
bool spacetime_plif::has_Ub() const {
	return my_Ub;
}

inline
const spacetime_plif::box_type& spacetime_plif::E_X0f(const duration& delta) const {
	return get_delta_params(delta).E_X0f;
}

inline
const spacetime_plif::box_type& spacetime_plif::E_X0b(const duration& delta) const {
	return get_delta_params(delta).E_X0b;
}

inline
const spacetime_plif::box_type& spacetime_plif::E_U(const duration& delta) const {
	return get_delta_params(delta).E_U;
}

inline
const support_function_provider::ptr& spacetime_plif::negAPhi2U(
		const duration& delta) const {
	return get_delta_params(delta).negAPhi2U;
}

inline
const spacetime_plif::box_type& spacetime_plif::E_Ub(const duration& delta) const {
	return get_delta_params(delta).E_Ub;
}

inline
const support_function_provider::ptr& spacetime_plif::negAPhi2Ub(
		const duration& delta) const {
	return get_delta_params(delta).negAPhi2Ub;
}

inline
const spacetime_plif::matrix_type& spacetime_plif::exp_AdeltaT(
		const duration& delta) const {
	return get_delta_params(delta).Phi_T;
}

inline
const spacetime_plif::matrix_type& spacetime_plif::phi1(const duration& delta) const {
	return get_delta_params(delta).Phi1;
}

inline
const spacetime_plif::delta_parameters& spacetime_plif::get_delta_params(
		const duration& delta) const {
	spacetime_plif* nonconst_this = const_cast<spacetime_plif*> (this);
	return nonconst_this->get_delta_params(delta);
}

inline
const spacetime_plif::delta_parameters& spacetime_plif::get_delta_params(
		const duration& delta) {
	if (min_infeasible_delta>=scalar_type(0) && delta>=min_infeasible_delta) {
		throw math::invalid_number_exception("throwing because of numerical issues, this should be caught");
	}

	if (my_delta_params_it == my_delta_params_cache.end() || !is_MEQ(
			my_delta_params_it->first, delta)) {
		my_delta_params_it = my_delta_params_cache.find(delta);
		if (my_delta_params_it == my_delta_params_cache.end()) {
			// not in the cache yet, so compute and add them
			try {
				delta_parameters params = compute_delta_params(delta);
//				std::cout << std::endl
//						<< "inserting delta=" << delta << " into cache" << std::endl;
				my_delta_params_it = my_delta_params_cache.insert_missing(delta,
						params);
//				std::cout << std::endl
//						<< "inserting done" << std::endl;

				min_delta = std::max(scalar_type(0),
						std::min(min_delta, delta));
				IFLOGGER(DEBUG5) {
					LOGGER(DEBUG5, "get_delta_params",
							"computed flowpipe parameters for delta="+to_string(delta));
				}
			} catch (math::invalid_number_exception& e) {
				if (min_infeasible_delta < scalar_type(0)) {
					min_infeasible_delta = delta;
				} else {
					min_infeasible_delta = std::min(min_infeasible_delta,
							delta);
				}
//				std::cout << std::endl
//						<< "throwing because of numerical issues" << std::endl;
				throw e;
			}

		}
	}
	return my_delta_params_it->second;
}

inline
spacetime_plif::delta_parameters spacetime_plif::compute_delta_params(
		const duration& delta) const {
	delta_parameters p;
	/* New computation of Phi at the same time as
	 * special matrices
	 * Phi=e^(A*delta)
	 * Phi1=A^{-1}(e^(A*delta)-I)
	 * Phi2=A^{-2}(e^(A*delta)-I-delta*A)
	 */
	const matrix_type& A = my_aff.dyn.get_A();
	positional_vdomain dom = A.domain();
	p.Phi = matrix_exponential(delta * A);
	p.Phi_T = p.Phi.transpose(); // Phi_T is the transpose of Phi.

	matrix_type A_abs = A.abs();
	matrix_type Phi_abs = matrix_type(dom, dom); // not used
	matrix_type Phi1_abs = matrix_type(dom, dom);
	matrix_type Phi2_abs = matrix_type(dom, dom);
	//matrix_type M2; // = A-A*exp(\delta*A)
	matrix_type M4; // = A^2*exp(\delta*A)
	matrix_type Mtemp;

	math::get_special_matrices(A_abs.get_matrix(), delta, Phi_abs.get_matrix(),
			Phi1_abs.get_matrix(), Phi2_abs.get_matrix());

	Mtemp = A * p.Phi;
	//M2 = A - Mtemp;
	M4 = A * Mtemp;

	vector_type zeros(dom);
	if (my_aff.X0) {
		support_function_provider const& X0 = *my_aff.X0;
		box_type BX2f = finite_bounding_box(my_S2f, to_map(Phi2_abs));
		//box_type S2b = finite_symmetric_bounding_box(X0, to_map(M4));
		box_type S2b = finite_symmetric_bounding_box(my_S2f, to_map(p.Phi));
		box_type BX2b = finite_bounding_box(S2b, to_map(Phi2_abs));
		p.E_X0f = BX2f;
		p.E_X0b = BX2b;
	} else {
		p.E_X0f = finite_hyperbox<scalar_type> (zeros, zeros);
		p.E_X0b = finite_hyperbox<scalar_type> (zeros, zeros);
	}

//	std::cout  << std::endl << "parameters for delta=" << delta << std::endl << " Phi2_abs:" << Phi2_abs <<std::endl <<" E_X0f:" << p.E_X0f <<std::endl << " E_X0b: " << p.E_X0b << std::endl;

	/** @todo Room for improvement:
	 * If possible, p.E_U and p.negAPhi2U could point to the same object
	 * as p.E_Ub and p.negAPhi2Ub, respectively.
	 */

	// we need Phi2 for the error bounds
	p.Phi1 = matrix_type(dom, dom); // not used
	p.Phi2 = matrix_type(dom, dom);
	if (has_Ub()) {
		math::get_special_matrices(A.get_matrix(), delta, p.Phi.get_matrix(),
				p.Phi1.get_matrix(), p.Phi2.get_matrix());
	}
	if (my_aff.dyn.get_U()) {
		box_type SU = finite_symmetric_bounding_box(*my_aff.dyn.get_U(), to_map(A));
		p.E_U = finite_bounding_box(SU, to_map(Phi2_abs));

		p.negAPhi2U = support_function_provider::ptr(
				new support_function::sf_unary<scalar_type>(my_aff.dyn.get_U(),
						to_map(-A * p.Phi2)));
	} else {
		// there is no U
		p.E_U = finite_hyperbox<scalar_type> (zeros, zeros);
		p.negAPhi2U = support_function_provider::ptr();
	}

	// now U+b
	if (my_Ub) {
		box_type SUb = finite_symmetric_bounding_box(*my_Ub, to_map(A));
		p.E_Ub = finite_bounding_box(SUb, to_map(Phi2_abs));
		p.negAPhi2Ub = support_function_provider::ptr(
				new support_function::sf_unary<scalar_type>(my_Ub,
						to_map(-A * p.Phi2)));
	} else {
		p.E_Ub = finite_hyperbox<scalar_type> (zeros, zeros);
		p.negAPhi2Ub = support_function_provider::ptr();
	}

	return p;
}

}

#endif /* SPACETIME_PLIF_DELTA_PARAMS_HPP_ */
