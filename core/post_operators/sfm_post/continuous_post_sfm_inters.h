/*
 * continuous_post_sfm_inters.h
 *
 *  Created on: Jan 20, 2010
 *      Author: frehse
 */

#ifndef CONTINUOUS_POST_SFM_INTERS_H_
#define CONTINUOUS_POST_SFM_INTERS_H_

#include "core/post_operators/continuous_post.h"

#include <typeinfo>

//#include "../abstract_framework/global_types.h"
//#include "../utility/shared_ptr_output.h"
//#include "../valuation_function/node_print_visitor.h"
#include <fstream>
#include "core/continuous/continuous_dynamics/ode_affine_dynamics.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/hybrid_automata/location.h"
#include "core/post_operators/sfm_post/support_function.h"

#include "core/continuous/support_function/sfm/sfm_cont_set.h"

namespace hybrid_automata {

/**
 * Computes the set of symbolic states that is reachable from cset with elapse of time upto delta.
 *
 * This version takes as directions only the negative invariant directions.
 *
 * \param csup The initial set from which the time elapse reachable set is computed
 * \param M The affine dynamics x'=Ax+b.
 * \param inv The invariant associated with the location.
 * \param time_step The time step.
 * \param time_horizon The time bound.
 *
 * Returns the set of symbolic states which contains the reachable time elapse set from cset. By symbolic state,
 * we mean a pair of discrete set and a continuous set. sstate_list is a collection of such symbolic states in which
 * the discrete set contains only one integer element k and the continuous set is a subset of R^n which is reachable
 * in the kth iteration of the discrete time elapse set computation algorithm. Here, the discrete time elapse computation
 * algorithm is the discrete version of the reachability algorithm using support function [Girard et al].
 */
template<class scalar_type> continuous::continuous_set::ptr support_function_time_elapse_inters(
		const continuous::support_function_provider::const_ptr& csup,
		const math::affine_map<scalar_type>& M,
		const typename continuous::polyhedron<scalar_type>::const_ptr& inv,
		double time_horizon, double time_step) {

	if (time_horizon < 0.0) {
		throw std::runtime_error(
				"support_function_time_elapse : cannot handle negative time_horizon "
						+ to_string(time_horizon));
	}
	if (time_step <= 0.0) {
		throw std::runtime_error(
				"support_function_time_elapse : cannot handle negative time_step "
						+ to_string(time_step));
	}

	if (csup->is_empty()) {
		continuous::continuous_set::ptr empty_set(csup->create_empty());
		return empty_set;
	}

	// call implementation of iterative computation of time elapse set using support function.

	//	typedef Rational scalar_type;
	typedef math::matrix<scalar_type> matrix_type;
	typedef math::vector<scalar_type> vector_type;

	//	dimension_t dim = csup->get_dim();


	unsigned int max_iterations;
	max_iterations = (unsigned int) ceil(time_horizon / time_step);

	std::list<vector_type> directions;
	// don't add any directions, just take the negative invariant
	// choose_directions(csup, directions); // R is the number of directions.

	matrix_type A = M.get_A().get_matrix();
	math::vdom_vector<scalar_type> b_named=M.get_b();

	// Create the set of inputs from b
	typename continuous::polyhedron<scalar_type>::ptr U =
			typename continuous::polyhedron<scalar_type>::ptr(new continuous::constr_polyhedron<
					scalar_type> ());
	add_vertice_constraints(*U, b_named);

	index_to_variable_id_map_ptr iimap = b_named.get_index_to_variable_id_map(); // iimap that tells which variables are at which index

	//	std::cout << "supp time elapse with dynamics A=" << A << " and U={" << U
	//			<< "}" << std::endl;
	matrix_type symbolic_reach_set;
	//cout << csup << endl;
	//cout << U << endl;

	// Add the inv directions to the directions_list. Also construct a vector of the inhomogenous terms of the invariants.
	// The order of the inhomogenous terms in the vector should be the same as the corresponding order of directions in the directions list.
	// For directions which are not invariant directions, the corresponding inhomogenous term is set to positive infinity.

	typename math::lin_constraint_system<scalar_type>::const_ptr
			inv_constr_ptr = inv->get_constraints();

	std::list<scalar_with_infinity<scalar_type> > inh_terms;

	scalar_with_infinity<scalar_type> my_scalar;

	/* This code below checks if a direction is one of the negative invariant directions. If it is not, then pos_infiity
	 * is stored for its corresponding entry in the inh_coeffs list. If the direction is one of the invariant
	 * direction, then the inhomogenous term of the invariant constraint is stored in the corresponding entry
	 * of inh_coeff list.
	 */
	for (typename std::list<vector_type>::const_iterator it =
			directions.begin(); it != directions.end(); it++) {
		bool flag = false;
		for (typename math::lin_constraint_system<scalar_type>::const_iterator
				iter = inv_constr_ptr->begin(); iter != inv_constr_ptr->end(); ++iter) {
			math::vdom_vector<scalar_type> rvec = iter->get_canonic_l();
			if (iimap != rvec.get_index_to_variable_id_map()) {
				rvec.reorder(iimap);
			}
			vector_type inv_direction = rvec.get_vector();

			//			std::cout << "inv_direction:" << inv_direction << std::endl;
			//			std::cout << "*it:" << *it << std::endl;

			if (*it == -inv_direction) {
				inh_terms.push_back(-scalar_with_infinity<scalar_type> (
						iter->get_canonic_inh_coeff()));
				flag = true; // The assumption is the inv set contains no redundant constraints
				break;
			}
		}
		if (!flag)
			inh_terms.push_back(my_scalar.pos_infty());
	}
	//Variant with Inv directions and not negative inv directions


	//	std::cout << "Invariant Set is:" << inv <<std::endl;

	// Variant with negative inv directions @author Rajarshi
	/*
	 * This code below checks if an negative invariant direction is already present in existing list of directions.
	 * Only the new negative invariant directions are added to the directions list.
	 */
	for (typename math::lin_constraint_system<scalar_type>::const_iterator
			it = inv_constr_ptr->begin(); it != inv_constr_ptr->end(); ++it) {
		bool red_flag = false;

		math::vdom_vector<scalar_type> rvec = it->get_canonic_l();
		if (iimap != it->get_l().get_index_to_variable_id_map()) {
			rvec.reorder(iimap);
		}
		vector_type inv_direction = rvec.get_vector();
		//		std::cout << "invariant direction is:" << inv_direction << std::endl;


		for (typename std::list<vector_type>::const_iterator iter =
				directions.begin(); iter != directions.end(); iter++) {
			if (*iter == -inv_direction) {
				red_flag = true; // the neg invariant is already present as a direction
				break;
			}
		}
		if (red_flag)
			continue;
		else {
			directions.push_back(-inv_direction); // Add the non-redundant negative inv direction to the list of directions
			inh_terms.push_back(-scalar_with_infinity<scalar_type> (
					it->get_canonic_inh_coeff()));
		}

	}

	//	testing
	/*	std::cout << "directions after adding inv directions : "<< std::endl;
	 for(typename std::list<vector_type>::const_iterator it = directions.begin(); it != directions.end(); it++){
	 std::cout << "direction:" << *it << std::endl;
	 }
	 */
	//make a vector from the list
	std::vector<scalar_with_infinity<scalar_type> >
			inh_coeffs(inh_terms.size());
	//	std::cout << "the inh coeffs are: " << std::endl;
		unsigned int i = 0;
		for (typename std::list<scalar_with_infinity<scalar_type> >::const_iterator
				it = inh_terms.begin(); it != inh_terms.end(); ++it, i++) {
			inh_coeffs[i] = *it;
			//		std::cout << *it << std::endl;
		}

	symbolic_reach_set
			= continuous::support_function::sf_evaluator<scalar_type>(A,
					b_named.get_vector(), iimap, csup, U).compute(directions,
					inh_coeffs, scalar_type(time_step), max_iterations);
	//cout << "Reach Set Matrix: " << symbolic_reach_set;

	// create the postc  parameter structure .
	continuous::support_function::postc_params<scalar_type> my_postc_params;
	my_postc_params.dynamics_A = A;
	my_postc_params.dynamics_b = b_named.get_vector();
	my_postc_params.initial_set_ptr = csup;
	my_postc_params.invariant_set_ptr = inv;
	my_postc_params.input_set_ptr = U;
	my_postc_params.delta = scalar_type(time_step);
	my_postc_params.time_horizon = scalar_type(time_horizon);

	continuous::support_function::sfm_cont_set<scalar_type>* sfm_ptr = new continuous::support_function::sfm_cont_set<scalar_type> (
			my_postc_params, symbolic_reach_set, directions, iimap);

	// @todo What is this for?
	//sfm_ptr->embed_variables(csup->get_variable_ids());
	//testing..
	//cout << "IIMAP :" << iimap << std::endl;

	return continuous::continuous_set::ptr(sfm_ptr);
}
;

/** \brief A continuous post operator for ode_affine_dynamics.
 *
 * This version is tailored for the precise intersection algorithm;
 * it only computes the SFM for the negative invariant directions.
 * The post_discrete operator is expected to request more directions
 * to the SFM created here.
 *
 * It is applicable under the following conditions:
 * - the dynamics are of type ode_affine_dynamics,
 * - the continuous_set class is a support_function_provider.
 */

template<class scalar_type>
class continuous_post_sfm_inters: public continuous_post {
public:
	typedef boost::shared_ptr<continuous_post_sfm_inters> ptr;
	typedef boost::shared_ptr<const continuous_post_sfm_inters> const_ptr;

	virtual ~continuous_post_sfm_inters() {
	}
	;
	virtual continuous_post_sfm_inters* clone() {
		return new continuous_post_sfm_inters(*this);
	}
	;

	/** Return the continuous_set that results from applying time elapse
	 * with the time constraints tcons to the continuous set *cset. */
	virtual continuous::continuous_set_ptr post(const time_constraints& tcons,
			const continuous::continuous_set_const_ptr& cset) const {
		typedef continuous::ode_affine_dynamics<scalar_type> ode_aff_dyn;

		continuous::continuous_set_ptr ret_set;
		if ( continuous::support_function_provider::const_ptr csup =
				boost::dynamic_pointer_cast<const continuous::support_function_provider > (cset)) {
			if (typename continuous::polyhedron<scalar_type>::const_ptr inv_poly_ptr =
				boost::dynamic_pointer_cast<
							const continuous::polyhedron<scalar_type> > (
							tcons.get_invariant())) {

				if (const ode_aff_dyn* dp=dynamic_cast<const ode_aff_dyn*>(tcons.get_dynamics().get())) {
					std::cout << "cont_post_inters "<< get_sampling_time() << std::endl;

					// call time elapse calculation
					// @todo There appears bugs with sf_time_elapse_inters, so change to
					// old time elapse.
//					ret_set = support_function_time_elapse_inters(csup, *dp,
//							inv_poly_ptr, get_time_horizon(),
//							get_sampling_time());
					ret_set = continuous::support_function::support_function_time_elapse(csup, *dp,
							inv_poly_ptr, get_time_horizon(),
							get_sampling_time());

					// Intersect the result with the invariant.
					//ret_set = compute_or_assign_intersection(ret_set,
					//		tcons.get_invariant());
										//debug
										typename continuous::support_function::sfm_cont_set<scalar_type>::ptr sfm_ptr =
												boost::dynamic_pointer_cast<continuous::support_function::sfm_cont_set<scalar_type> >(ret_set);

										std::cout << "sfm after inv intersection" << std::endl;
										continuous::polyhedron<double>::set_output_format(continuous::polyhedron<double>::DOUBLE_GENERATORS);
										std::ofstream myfile;
										myfile.open ("/tmp/out_tmp.tmp");

										for(unsigned int i=0;i<sfm_ptr->get_sfm().size2();i++)
										{
											continuous::constr_polyhedron<double> poly = sfm_ptr->get_polytope(i);
											myfile << poly; myfile << std::endl << std::endl;
										}
										myfile.close();
										system("graph -TX -C -B -q0.5 /tmp/out_tmp.tmp");
										//end debug


				} else {
					//				std::cerr << "offended by dynamics:"
					//						<< tcons.get_dynamics()->get_predicate() << std::endl;
					std::string tname =
							typeid(*tcons.get_dynamics().get()).name();
					throw std::runtime_error(
							"continuous_post_sfm_inters::post cannot handle dynamics of type "
									+ tname);
				}
			} else {
				std::string tname = typeid(*tcons.get_invariant().get()).name();
				throw std::runtime_error(
						"continuous_post_sfm_inters::post cannot handle invariant of type "
								+ tname);
			}
		} else {
			std::string tname = typeid(*cset.get()).name();
			throw std::runtime_error(
					"continuous_post_sfm_inters::post cannot handle continuous_set of type"
							+ tname);
		}

		return ret_set;
	}
	;

};

}

#endif /* continuous_post_sfm_INTERS_H_ */
