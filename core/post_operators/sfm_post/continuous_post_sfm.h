/*
 * continuous_post_sfm.h
 *
 *  Created on: Nov 12, 2009
 *      Author: frehse
 */

#ifndef CONTINUOUS_POST_SFM_H_
#define CONTINUOUS_POST_SFM_H_

#include "utility/logger_stopwatch.h"

#include "core/post_operators/continuous_post.h"

#include <typeinfo>

//#include "../abstract_framework/global_types.h"
//#include "../utility/shared_ptr_output.h"
//#include "../valuation_function/node_print_visitor.h"
#include "core/continuous/continuous_dynamics/ode_affine_dynamics.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/hybrid_automata/location.h"
#include "core/post_operators/sfm_post/support_function.h"
#include "core/continuous/support_function/sf_base/sf_unary.h"

namespace hybrid_automata {

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
class continuous_post_sfm: public continuous_post {
public:
	typedef boost::shared_ptr<continuous_post_sfm> ptr;
	typedef boost::shared_ptr<const continuous_post_sfm> const_ptr;

	continuous_post_sfm(double flowpipe_tolerance) :
		my_flowpipe_tolerance(flowpipe_tolerance) {
	}
	;

	virtual ~continuous_post_sfm() {
	}
	;

	virtual continuous_post_sfm* clone() {
		return new continuous_post_sfm(*this);
	}
	;

	/** Return the continuous_set that results from applying time elapse
	 * with the time constraints tcons to the continuous set *cset. */
	virtual continuous::continuous_set_ptr post(const time_constraints& tcons,
			const continuous::continuous_set_const_ptr& cset) const {
		LOGGERSW(DEBUG1,"continuous_post_sfm::post","Continuous post with continuous_post_sfm");

		using namespace continuous;

		typedef ode_affine_dynamics<scalar_type> ode_aff_dyn;

		continuous_set_ptr ret_set = continuous_set_ptr();
		if ( support_function_provider::const_ptr csup =
				boost::dynamic_pointer_cast<const support_function_provider > (cset)) {
			typename polyhedron<scalar_type>::const_ptr inv_poly_ptr;
			if (tcons.get_invariant()) {
				inv_poly_ptr = boost::dynamic_pointer_cast<const polyhedron<
						scalar_type> >(tcons.get_invariant());
				if (!inv_poly_ptr) {
					std::string tname =
							typeid(*tcons.get_invariant().get()).name();
					throw basic_exception(
							"continuous_post_sfm::post cannot handle invariant of type "
									+ tname);
				}
			}

			if (const ode_aff_dyn* dp=dynamic_cast<const ode_aff_dyn*>(tcons.get_dynamics().get())) {
				LOGGER(DEBUG3,"continuous_post_sfm::post","computing support_function_time_elapse");

				// first, we treat degenerated cases
				if (dp->is_universe()) {
					// there is no restriction on the dynamics,
					// so the entire invariant is reachable.
					// simply return a copy of the invariant

					LOGGER(DEBUG3,"continuous_post_sfm::post","no flow restriction (flow predicate is always true)");
					ret_set = continuous_set_ptr(inv_poly_ptr->clone());
				} else if (dp->is_zero()) {
					// there is no change in the variables,
					// since the derivatives are identical zero.
					// Simply return a copy of the initial states

					LOGGER(DEBUG3,"continuous_post_sfm::post","no time elapse in location (unsat flow predicate)");
					ret_set = continuous_set_ptr(cset->clone());
				} else {
					// call time elapse calculation
					// restrict invariant with const variables
					typename polyhedron<scalar_type>::const_ptr new_inv =
							inv_poly_ptr;
					if (false && new_inv) {
						const math::affine_map<scalar_type>& M_orig(*dp);
						// First, construct the bounding box of the derivatives in this invariant
						support_function::sf_unary<scalar_type> deriv_set(
								new_inv, M_orig);
						// Go through the variables in M_orig and find the ones with derivative zero
						variable_id_set constvars = get_vars_bound_to_zero<
								scalar_type> (deriv_set);

						if (!constvars.empty()) {
							std::stringstream ss;
							logger::copyfmt_to(ss);
							print_variable_id_set(ss, constvars);
							LOGGER(DEBUG4,"continuous_post_sfm","adding const constraints to variables "+ss.str());
							//std::cout << "constant state vars:"; print_variable_id_set(std::cout,constvars); std::cout << std::endl;
							// Get bounds on the initial values of constvars
							positional_vdomain const_dom(constvars);
							typedef math::vector<scalar_type> vector_type;
							std::list<vector_type> constdirections;
							support_function::choose_directions(const_dom, constdirections); // R is the number of directions.
							typedef std::set<math::vdom_vector<scalar_type>,
									math::numeric::lex_comp_less<scalar_type,
											math::vdom_vector> > dir_set_type;

							dir_set_type dir_set;

							for (typename std::list<math::vector<scalar_type> >::const_iterator
									it = constdirections.begin(); it
									!= constdirections.end(); ++it) {
								math::vdom_vector<scalar_type>
										d(const_dom, *it);
								dir_set.insert(d);
							}

							typename constr_polyhedron<scalar_type>::ptr poly(
									new constr_polyhedron<scalar_type> (
											compute_outer_poly<scalar_type> (
													*csup, dir_set)));
							IFLOGGER(DEBUG5) {
								ss.str("");
								ss << poly;
								LOGGER(DEBUG5,"continuous_post_sfm","const constraints are: "+ss.str());
							}
							// add poly to new_inv
							poly->add_constraints(new_inv->get_constraints());
							new_inv = poly;
							//std::cout << "new invariant:" << new_inv;
						}
					}

					typename support_function::sfm_cont_set<scalar_type>::ptr
							sfm_ptr;
					support_function_provider::const_ptr U = boost::static_pointer_cast<
							const support_function_provider>(dp->get_U());
					sfm_ptr = support_function::support_function_time_elapse(
							csup, *dp, U, new_inv, get_time_horizon(),
							get_sampling_time(), my_flowpipe_tolerance);
					LOGGER(DEBUG3,"continuous_post_sfm::post","intersecting flowpipe ("+to_string(sfm_ptr->get_size())+" steps, "+to_string(sfm_ptr->get_directions().size())+" directions) with invariant");

					// Intersect the result with the invariant
					// The result is not empty because the root of the sfm is not empty.
					if (new_inv) {
						sfm_ptr->intersection_with_poly(*new_inv);
					}
					ret_set = sfm_ptr;
				}
			} else {
				//				std::cerr << "offended by dynamics:"
				//						<< tcons.get_dynamics()->get_predicate() << std::endl;
				std::string tname = typeid(*tcons.get_dynamics().get()).name();
				throw std::runtime_error(
						"continuous_post_sfm::post cannot handle dynamics of type "
								+ tname);
			}
		} else {
			std::string tname = typeid(*cset.get()).name();
			throw std::runtime_error(
					"continuous_post_sfm::post cannot handle continuous_set of type"
							+ tname);
		}

		return ret_set;
	}
	;

private:
	double my_flowpipe_tolerance;
};

}

#endif /* CONTINUOUS_POST_SFM_H_ */
