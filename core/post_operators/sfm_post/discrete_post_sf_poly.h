/*

 * discrete_post_sf_poly.h
 *
 *  Created on: Oct 17, 2010
 *      Author: frehse
 */

#ifndef DISCRETE_POST_SF_POLY_H_
#define DISCRETE_POST_SF_POLY_H_

#include <typeinfo>

#include "utility/logger_stopwatch.h"

#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "core/hybrid_automata/transition.h"
#include "core/post_operators/discrete_post.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/continuous/support_function/sf_derived/sf_sum.h"
#include "core/continuous/support_function/sf_base/sf_unary_utility.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/polyhedra/polyhedron_collection.h"
#include "core/continuous/polyhedra/polyhedron_collection_operators.h"
#include "core/continuous/polyhedra/constraint_hull_generator.h"

namespace hybrid_automata {

/** A discrete post operator for sfm.
 *
 * The intersection computation is applied to the template hull of
 * the intersecting polyhedra. */
template<class scalar_type>
class discrete_post_sf_poly: public discrete_post {
public:
	typedef boost::shared_ptr<discrete_post_sf_poly<scalar_type> > ptr;
	typedef boost::shared_ptr<const discrete_post_sf_poly<scalar_type> >
			const_ptr;

	typedef typename continuous::support_function::sfm_cont_set<scalar_type>
			sfm;
	typedef typename sfm::ptr sfm_ptr;
	typedef typename sfm::const_ptr sfm_const_ptr;
	typedef typename continuous::polyhedron<scalar_type>::const_ptr
			poly_const_ptr;

	typedef continuous::polyhedron_collection<scalar_type> poly_list_type;
	typedef typename poly_list_type::const_ptr poly_list_const_ptr;

	discrete_post_sf_poly(double aggregation_type, double clustering) :
		my_aggregation_type(aggregation_type), my_clustering(clustering) {
	}
	;

	virtual ~discrete_post_sf_poly() {
	}
	;

	discrete_post_sf_poly* clone() {
		return new discrete_post_sf_poly<scalar_type> (*this);
	}
	;

	struct problem_instance {
		sfm_const_ptr sfm_set_ptr;
		poly_list_const_ptr cset_poly_coll_ptr;
		poly_const_ptr cset_poly_ptr;
		poly_const_ptr sinv_ptr;
		poly_const_ptr tinv_ptr;
		poly_const_ptr guard_poly_ptr;
		const typename continuous::reset_affine_transform<scalar_type>
				* reset_trans_ptr;
	};

	problem_instance get_instance(const jump_constraints& j,
			continuous::continuous_set::const_ptr source_inv,
			continuous::continuous_set::const_ptr target_inv,
			continuous::continuous_set::const_ptr cset) const {
		problem_instance inst;
		inst.sfm_set_ptr
				= boost::dynamic_pointer_cast<
						const continuous::support_function::sfm_cont_set<
								scalar_type> >(cset);
		if (!inst.sfm_set_ptr) {
			inst.cset_poly_coll_ptr = boost::dynamic_pointer_cast<
					const poly_list_type>(cset);
			if (!inst.cset_poly_coll_ptr) {
				inst.cset_poly_ptr = boost::dynamic_pointer_cast<
						const continuous::polyhedron<scalar_type> >(cset);
				if (!inst.cset_poly_ptr) {
					std::string tname = typeid(*cset.get()).name();
					throw std::runtime_error(
							"discrete_post_sf_poly: continuous_set type "
									+ tname + " not supported as states");
				}
			}
		}
		inst.sinv_ptr = boost::dynamic_pointer_cast<
				const continuous::polyhedron<scalar_type> >(source_inv);
		if (source_inv && !inst.sinv_ptr) {
			throw std::runtime_error(
					"discrete_post_sf_poly: continuous_set type not supported as invariant");
		}
		inst.tinv_ptr = boost::dynamic_pointer_cast<
				const continuous::polyhedron<scalar_type> >(target_inv);
		if (target_inv && !inst.tinv_ptr) {
			throw std::runtime_error(
					"discrete_post_sf_poly: continuous_set type not supported as invariant");
		}
		inst.guard_poly_ptr = poly_const_ptr();
		if (j.get_guard()) {
			inst.guard_poly_ptr = boost::dynamic_pointer_cast<
					const continuous::polyhedron<scalar_type> >(j.get_guard());
			if (!inst.guard_poly_ptr) {
				throw std::runtime_error(
						"discrete_post_sf_poly: guard type not supported");
			}
		}
		inst.reset_trans_ptr = 0;
		if (j.get_transform()) {
			inst.reset_trans_ptr
					= dynamic_cast<const continuous::reset_affine_transform<
							scalar_type>*> (j.get_transform().get());
			if (!inst.reset_trans_ptr) {
				std::stringstream ss;
				logger::copyfmt_to(ss);
				ss << "offending transform:" << j.get_transform();
				throw std::runtime_error(
						"discrete_post_sf_poly: transform type not supported\n"
								+ ss.str());
			}
		}
		return inst;
	}
	;

	void post_poly(problem_instance& inst, poly_list_type & polys,
			variable_id_set polys_vis) const {

		LOGGERSWOC(DEBUG4,__FUNCTION__,"Computing discrete successors of polyhedra");

		using namespace continuous;
		using namespace support_function;
		if (inst.reset_trans_ptr) {
			// a pointer for the inputs
			support_function_provider::const_ptr mapped_U =
					boost::static_pointer_cast<const support_function_provider>(
							inst.reset_trans_ptr->get_input_set());

			math::affine_map<scalar_type> M = *inst.reset_trans_ptr;

			// Remove variables not in the domain of M
			// get the variables from invariant, guard and sfm
			variable_id_set excess_variables = polys_vis;
			set_difference_assign(excess_variables,
					M.domain().get_variable_ids());
			if (!excess_variables.empty()) {
				if (quantify_unused_jump_variables) {
					// Use existential quantification to eliminate variables
					LOGGER(DEBUG5,"discrete_post_f_poly::post","quantifying over "+to_string(excess_variables.size())+" variable(s) not in assignment");

					// additional log output
					IFLOGGER(DEBUG7) {
						std::stringstream ss;
						logger::copyfmt_to(ss);
						print_variable_id_set(ss, polys_vis);
						LOGGER_OS(DEBUG7,"discrete_post_f_poly::post")
							<< "poly variables: " << ss.str() << ", ";
						ss.str("");
						print_variable_id_set(ss, M.domain().get_variable_ids());
						LOGGER_OS(DEBUG7,"discrete_post_f_poly::post")
							<< "map domain: " << ss.str() << ", " << "affine map: "
									<< M << std::flush;
					}

					bool bef_empty = polys.is_empty();
					//std::cout << "before quant:"<< bef_empty << ", ";
					poly_list_type qpolys(*polys.clone());
					qpolys.existentially_quantify_variables(excess_variables);
					bool aft_empty = qpolys.is_empty();
					//std::cout << "after quant:"<< aft_empty;
					if (false && bef_empty != aft_empty) {
						std::stringstream ss;
						logger::copyfmt_to(ss);
						ss << "before: " << polys << std::endl;
						ss << "after: " << qpolys;
						constr_polyhedron<scalar_type> testpoly;
						testpoly.add_constraints(
								*(*polys.begin())->get_constraints());
						ss << "test: " << testpoly.is_empty();
						throw basic_exception(
								"internal error: emptiness changed due to quantification.\n"
										+ ss.str());
					}
					polys = qpolys;
				} else {
					// add the variables to the domain of the map

					positional_vdomain new_dom=M.domain();
					new_dom.add_variables(create_variable_set(excess_variables));
					M.reorder(M.codomain(),new_dom);

					IFLOGGER(DEBUG7) {
						// additional log output
						LOGGER_OS(DEBUG7,"discrete_post_f_poly::post") << "Changed domain of jump map to include variables of continuous set: ";
						std::stringstream ss;
						logger::copyfmt_to(ss);
						print_variable_id_set(ss, polys_vis);
						LOGGER_OS(DEBUG7,"discrete_post_f_poly::post") << ss.str() << ". New matrix:" << M << std::endl << "new codomain:"<<M.codomain()<<", new domain:" << M.domain() << std::endl;
					}
				}
			}

			// If the domain and codomain are not the same,
			// we need to reorder and seperate state from input variables
			if (M.domain() != M.codomain()) {
				// split A into square A1 and A2
				math::affine_map<scalar_type> Mstate;
				math::affine_map<scalar_type> Minput;
				math::separate_states_from_inputs(M, Mstate, Minput);

				M = Mstate;
//									std::cout << "split dynamics into " << Mstate << " and " << Minput
//											<< " with input domain " << Minput.domain() << std::endl;

				// Define the set of inputs using the invariant.
				// This is the original invariant with both state and nonstate variables
				if (Minput.domain().size() > 0) {
					LOGGER(DEBUG5,"discrete_post_f_poly::post","nondeterministic inputs "+logger::formatted_to_string(Minput.domain()));
					if (inst.reset_trans_ptr->get_input_set()) {
						LOGGER(DEBUG5,"discrete_post_f_poly::post","using input set "+logger::formatted_to_string(inst.reset_trans_ptr->get_input_set()));
						mapped_U = boost::static_pointer_cast<const support_function_provider>(inst.reset_trans_ptr->get_input_set());
					} else if (inst.sinv_ptr) {
						LOGGER(DEBUG5,"discrete_post_f_poly::post","using source invariant "+logger::formatted_to_string(inst.sinv_ptr));
						mapped_U = inst.sinv_ptr;
					} else {
						// make a universe set for this domain
						mapped_U = support_function_provider::const_ptr(
								new hyperbox<scalar_type> ());
					}
					mapped_U = support_function_provider::const_ptr(
							new support_function::sf_unary<scalar_type>(
									mapped_U, Minput));
				}
			}

			// If A is singular or there are nondet. inputs,
			// map the support function and the obtain an outer poly approx.
			// Otherwise, map the polys directly (exact solution).
			scalar_type det = math::matrix_determinant(M.get_A().get_matrix());
			bool singular = math::numeric::is_MEQ(det / scalar_type(1000),
					scalar_type(0));
			if (singular || mapped_U) {
				if (singular && !mapped_U) {
					LOGGER(DEBUG4,"discrete_post_f_poly::post","using outer poly for mapping because of non-regular jump");
				} else if (!singular && mapped_U) {
					LOGGER(DEBUG4,"discrete_post_f_poly::post","using outer poly for mapping because of nondet. inputs");
				} else {
					LOGGER(DEBUG4,"discrete_post_f_poly::post","using outer poly for mapping because of non-regular jump and nondet. inputs");
				}
				// ... use sf to transform and get outer poly

				// get the directions for the outer polys
				std::list<math::vector<scalar_type> > directions;
				choose_directions(M.codomain(), directions); // R is the number of directions.
				typename sf_set<scalar_type>::vector_set sf_directions;
				for (typename std::list<math::vector<scalar_type> >::const_iterator
						it = directions.begin(); it != directions.end(); ++it) {
					typename sf_set<scalar_type>::vector_type sf_vec(
							M.codomain(), *it);
					sf_directions.insert(sf_vec);
				}

				// make a new list of polys with the result
				poly_list_type mapped_polys;
				for (typename poly_list_type::const_iterator it = polys.begin(); it
						!= polys.end(); ++it) {
					// the result is the minkowski sum of the mapped states and the mapped inputs

					// get the mapped states
					typename sf_set<scalar_type>::ptr
							mapped_states =
									typename sf_set<scalar_type>::ptr(
											new support_function::sf_unary<
													scalar_type>(*it, M));

					// form the minkowski sum with the inputs
					typename sf_set<scalar_type>::ptr mapped_sf;
					if (mapped_U) {
						mapped_sf = typename sf_set<scalar_type>::ptr(
								new sf_sum<scalar_type> (mapped_states,
										mapped_U));
					} else {
						mapped_sf = mapped_states;
					}

					// get the outer poly of the mapped sf object
					typename sf_set<scalar_type>::poly_ptr mapped_poly =
							mapped_sf->outer_poly(sf_directions);
					// add the poly to the result
					mapped_polys.insert(mapped_poly);
				}

				// replace polys with the mapped polys
				polys.swap(mapped_polys);
			} else {
				LOGGERSW(DEBUG4,"discrete_post_f_poly::post","applying affine map to poly in discrete_post_f_poly");
				// transform the polys
				const math::affine_map<scalar_type>& M = *inst.reset_trans_ptr;
				polys = apply_map(polys, M);
			}
		}
		if (inst.tinv_ptr) {
			// add constraints with redundancy check because
			// non-redundancy is needed later for template hull
			//polys.add_constraints(*inst.tinv_ptr->get_constraints(), true);
			IFLOGGER(DEBUG5) {
				std::stringstream ss;
				logger::copyfmt_to(ss);
				ss << inst.tinv_ptr;
				LOGGER(DEBUG5,"discrete_post_f_poly::post","adding target invariant "+ss.str());
			}
			polys.add_constraints(*inst.tinv_ptr->get_constraints(), true);
			polys.remove_empty();
		}
		polys.remove_empty();
	}

	/** Current implementation works only for continuous sets as polytopes represented as SFM.
	 */
	continuous::continuous_set_collection post(const jump_constraints& j,
			continuous::continuous_set::const_ptr source_inv,
			continuous::continuous_set::const_ptr target_inv,
			continuous::continuous_set::const_ptr cset) const {
		using namespace continuous;
		using namespace support_function;
		LOGGERSW(DEBUG1,"discrete_post_f_poly::post","Discrete post");

		typedef math::numeric::interval<unsigned int> intv_type;
		typedef std::list<intv_type> intv_list_type;
		typedef constr_polyhedron<scalar_type> poly_type;
		typedef polyhedron_collection<scalar_type> poly_list_type;

		problem_instance inst = get_instance(j, source_inv, target_inv, cset);
		continuous_set_collection ret_collection;

		// add the pre-transformed invariant constraints to the guard constraint
		typename constr_polyhedron<scalar_type>::ptr guard_poly_ptr = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
		constr_polyhedron<scalar_type>& guard_poly(*guard_poly_ptr);
		if (inst.guard_poly_ptr){
			// add guard constraints (duh)
			guard_poly.add_constraints(inst.guard_poly_ptr->get_constraints());
		}
		if(inst.sinv_ptr){
			// add source invariant constraints
			guard_poly.add_constraints(inst.sinv_ptr->get_constraints());
		}
		if(inst.tinv_ptr){
			constr_polyhedron<scalar_type> exit_set_poly;
			typename constr_polyhedron<scalar_type>::const_ptr tgt_inv_poly_ptr =
					boost::dynamic_pointer_cast<const constr_polyhedron<scalar_type> >(inst.tinv_ptr);
			if(tgt_inv_poly_ptr){
				math::affine_map<scalar_type> M = *inst.reset_trans_ptr;

				support_function_provider::const_ptr U =
						boost::dynamic_pointer_cast<const support_function_provider>(
								inst.reset_trans_ptr->get_input_set());

				constr_polyhedron<scalar_type> dummy;
				reverse_map_constraints(tgt_inv_poly_ptr, M , U, exit_set_poly, dummy);
				LOGGER_OS(DEBUG5,"discrete_post_f_poly::post") << "invariant poly reverse mapped to: " << exit_set_poly;
				guard_poly.add_constraints(exit_set_poly.get_constraints());
			}
			else
				throw std::runtime_error("discrete_post_f_poly: Target location invariant set type not supported\n");
		}
		LOGGER_OS(DEBUG5,"discrete_post_f_poly::post") << "guard cons after rev map:" << guard_poly;
		guard_poly.remove_redundant_constraints();
		LOGGER_OS(DEBUG5,"discrete_post_f_poly::post") << "guard constraints after redundancy check: " << guard_poly;
		//collapse guard constraints
		guard_poly.collapse_inequalities();
		LOGGER_OS(DEBUG5,"discrete_post_f_poly::post") << "guard constraints list after collapsing : " << guard_poly;
		inst.guard_poly_ptr = guard_poly_ptr;

		if (!inst.reset_trans_ptr->is_void()) {

			// If the set is an sfm
			if (inst.sfm_set_ptr) {
				LOGGERSW(DEBUG2,"discrete_post_f_poly::post","computing discrete post of sfm");
				LOGGER_ATTACH(DEBUG2,"discrete_post_f_poly::post"," of size "+to_string(inst.sfm_set_ptr->get_size()),logger::get_last_id());

				// Initialize the list with all entries of the sfm
				intv_list_type list;
				if (inst.sfm_set_ptr->get_size() > 0) {
					intv_type I; // universe
					I.set_lower(0);
					I.set_upper(inst.sfm_set_ptr->get_size() - 1);
					list.push_back(I);
				}

				// Keep track of the variables in the polys
				variable_id_set polys_vis =
						inst.sfm_set_ptr->get_variable_ids();
				set_union_assign(polys_vis, inst.sinv_ptr->get_variable_ids());

				//std::cout << "initially list: " << list << std::endl;

				// @attention We don't intersect with the source invariant because we assume
				// it has already been dealt with.
				// This might change in the future!

				// Find intersecting intervals
				//sfm_const_ptr sfm2;
				if (inst.guard_poly_ptr) {
					LOGGERSWOC(DEBUG5,__FUNCTION__,"Getting guard interval");
					// need to const_cast because the sfm will be extended
					typename sfm_cont_set<scalar_type>::ptr
							nonconst_sfm = boost::const_pointer_cast<
									sfm_cont_set<scalar_type> >(
									inst.sfm_set_ptr);
					list = get_intersection_interval(nonconst_sfm,
							*inst.guard_poly_ptr);
					//std::cout << "getting list from guard: " << list << std::endl;
				}
				else
					LOGGER(DEBUG2,"discrete_post_f_poly::post","There is no guard in transition");

				LOGGER(DEBUG2,"discrete_post_f_poly::post","found "+to_string(list.size())+" intervals intersecting with guard");

				for (intv_list_type::const_iterator it = list.begin(); it
						!= list.end(); ++it) {
					//std::cout << "treating interval: " << *it << std::endl;
					LOGGER(DEBUG3,"discrete_post_f_poly::post","treating interval of size "+to_string(it->size()));


					// create a new poly collection
					typename poly_list_type::ptr polys_ptr =
							typename poly_list_type::ptr(new poly_list_type());
					poly_list_type& polys = *polys_ptr;

					{
						LOGGERSWOC(DEBUG5,__FUNCTION__,"Intersection with guard");

						// get outer polys for the current interval
						for (unsigned int i = it->lower().get_val(); i
								<= it->upper().get_val(); ++i) {
							// insert without redundancy check
							typename polyhedron<scalar_type>::ptr
									new_poly(
											new constr_polyhedron<scalar_type> (
													inst.sfm_set_ptr->get_polytope(
															i)));
//							new_poly->remove_redundant_constraints();
							polys.insert(new_poly, false);
						}

						// intersect with the guard
						// @todo this might be no longer necessary if the
						// directions of the guard were added to the sfm
						// while getting the intervals
						if (inst.guard_poly_ptr) {
							// add constraints and perform redundancy removal because
							// non-redundancy is needed later for template hull
							polys.add_constraints(
									*inst.guard_poly_ptr->get_constraints(),
									true);

							// update the variables
							set_union_assign(polys_vis,
									inst.guard_poly_ptr->get_variable_ids());
						}
						polys.remove_empty();
					}

					if (!math::definitely(polys.is_empty())) {
						post_poly(inst, polys, polys_vis);
					}

					polys.remove_redundant_constraints();

					// either pass on the convex hull of the polys
					// or the template hull...
					if (polys.size() > 1 && my_clustering > 0.0
							&& my_aggregation_type != 1) {
						double err = my_clustering;
						// there's no need to do this if aggregation is thull without size limit
						LOGGERSW(DEBUG3,"discrete_post_f_poly::post","clustering");
						LOGGER_ATTACH(DEBUG3,"discrete_post_f_poly::post"," "+to_string(polys.size()) + " polyhedra using template hull with factor "+to_string(err*100),logger::get_last_id());
						polys
								= constraint_hull_generator<scalar_type> (
										*polys_ptr).compute_constraint_hull_up_to_error(err);

						LOGGER_ATTACH(DEBUG3,"discrete_post_f_poly::post"," obtained "+to_string(polys.size()) + " polyhedral clusters, ",logger::get_last_id());
					}
					if (!math::definitely(polys.is_empty())) {
						if (my_aggregation_type == 0 || double(polys.size())
								<= my_aggregation_type) {
							// use convex hull, i.e., insert the poly_collection to be used as sf_provider

							// if it's just one set, unwrap
							if (polys.size() == 1) {
								LOGGER(DEBUG3,"discrete_post_f_poly::post","adding " + to_string(polys.size()) + " polyhedron to result");
								ret_collection.insert(*polys.begin());
							} else {
								LOGGER(DEBUG3,"discrete_post_f_poly::post","adding convex hull of " + to_string(polys.size()) + " polyhedra to result");
								//continuous::support_function_provider::ptr
								continuous_set::ptr
										chull_ptr =
												support_function::construct_convex_hull<
														scalar_type,
														poly_list_type>(
														polys_ptr);
								ret_collection.insert(chull_ptr);
							}
						} else if (my_aggregation_type > 0) {
							// use template hull
							LOGGERSW(DEBUG3,"discrete_post_f_poly::post","adding polyhedral template hull");
							LOGGER_ATTACH(DEBUG3,"discrete_post_f_poly::post"," of "+to_string(polys.size()) + " polyhedra to result, ",logger::get_last_id());
							ret_collection.insert(
									constraint_hull_generator<scalar_type> (
											*polys_ptr).compute_constraint_hull());
						} else {
							// treat each convex set separately
							//  -- creating a new symbolic state for each one
							LOGGER(DEBUG3,"discrete_post_f_poly::post","adding " + to_string(polys.size()) + " separate polyhedra to result");
							for (typename polyhedron_collection<scalar_type>::const_iterator
									it = polys.begin(); it != polys.end(); ++it) {
								ret_collection.insert(*it);
							}
						}
					} else {
						LOGGER(DEBUG3,"discrete_post_f_poly::post","result of jump is empty");
					}
				}
				//					std::cout << "guard intersect : " << guard_intersect << std::endl;
				//					std::cout << "transform : " << *reset_trans_ptr << std::endl;
			} else if (inst.cset_poly_coll_ptr || inst.cset_poly_ptr) {
				// compute the discrete post of a polyhedron

				// Keep track of the variables in the polys
				variable_id_set polys_vis = cset->get_variable_ids();
				set_union_assign(polys_vis, inst.sinv_ptr->get_variable_ids());

				//std::cout << "initially list: " << list << std::endl;

				// @attention We don't intersect with the source invariant because we assume
				// it has already been dealt with.
				// This might change in the future!

				LOGGERSW(DEBUG2,"discrete_post_f_poly::post","computing discrete post of polyhedron");

				// create a new poly collection
				// and add the poly to it
				poly_list_type polys;
				if (inst.cset_poly_coll_ptr) {
					// clone the polys
					for (typename polyhedron_collection<scalar_type>::const_iterator
							it = inst.cset_poly_coll_ptr->begin(); it != inst.cset_poly_coll_ptr->end(); ++it) {
						polys.insert(typename polyhedron<scalar_type>::ptr(
								(*it)->clone()));
					}
				} else {
					polys.insert(typename polyhedron<scalar_type>::ptr(
							inst.cset_poly_ptr->clone()));
				}

				// intersect with the guard
				// @todo this might be no longer necessary if the
				// directions of the guard were added to the sfm
				// while getting the intervals
				if (inst.guard_poly_ptr) {
					// add constraints with redundancy check because
					// non-redundancy is needed later for template hull
					polys.add_constraints(
							*inst.guard_poly_ptr->get_constraints(), true);
					//polys.add_constraints(*inst.guard_poly_ptr->get_constraints(),false);
					set_union_assign(polys_vis,
							inst.guard_poly_ptr->get_variable_ids());
				}
				polys.remove_empty();

				post_poly(inst, polys, polys_vis);

				for (typename polyhedron_collection<scalar_type>::const_iterator
						it = polys.begin(); it != polys.end(); ++it) {
					ret_collection.insert(*it);
				}
			}
		} else {
			// null transform pointer
			throw basic_exception("internal error: unexpected null pointer for discrete transform");
		}

		//		else {
		//				// the transform is empty
		//				continuous::continuous_set_ptr ret_set
		//				= continuous::continuous_set_ptr(
		//						new continuous::constr_polyhedron<scalar_type>(
		//								continuous::constr_polyhedron<scalar_type>::empty_poly()));
		//				ret_collection.insert(ret_set);
		//			}

		return ret_collection;
	}
	;

private:
	double my_aggregation_type;
	double my_clustering;

public:
	static bool quantify_unused_jump_variables;
};

template<class scalar_type>
bool discrete_post_sf_poly<scalar_type>::quantify_unused_jump_variables = true;

}

#endif /* DISCRETE_POST_SF_POLY_H_ */
