/*
 * discrete_post_intv_split.h
 *
 *  Created on: Sep 8, 2011
 *      Author: ray
 */

#ifndef DISCRETE_POST_INTV_SPLIT_H_
#define DISCRETE_POST_INTV_SPLIT_H_


#include "utility/logger_stopwatch.h"

#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "core/hybrid_automata/transition.h"
#include "core/post_operators/discrete_post.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/continuous/support_function/sf_derived/sf_sum.h"
#include "core/continuous/support_function/sf_derived/sf_chull.h"
#include "core/continuous/support_function/sfm/sfm_section.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/polyhedra/polyhedron_collection.h"
#include "core/continuous/polyhedra/polyhedron_collection_operators.h"
#include "core/continuous/polyhedra/constraint_hull_generator.h"
#include "core/continuous/support_function_provider_utility.h"
#include "core/continuous/support_function/sf_derived/sf_lin_con_intersection.h"
#include "math/vdom/positional_vdomain.h"
#include "io/GEN_format/GEN_formatter.h"
#include "core/continuous/continuous_set_operator_implementations/intersection.h"
#include "core/continuous/support_function/sfm/sfm_lin_con_intersection.h"

namespace hybrid_automata {

/** A discrete post operator for sfm.
 *
 */

template<class scalar_type>
class discrete_post_intv_split: public discrete_post {
public:
	typedef boost::shared_ptr<discrete_post_intv_split<scalar_type> > ptr;
	typedef boost::shared_ptr<const discrete_post_intv_split<scalar_type> >
			const_ptr;

	typedef typename continuous::support_function::sfm_cont_set<scalar_type>
			sfm;
	typedef typename sfm::ptr sfm_ptr;
	typedef typename sfm::const_ptr sfm_const_ptr;
	typedef typename continuous::polyhedron<scalar_type>::const_ptr
			poly_const_ptr;
	typedef continuous::polyhedron_collection<scalar_type> poly_list_type;
	typedef typename poly_list_type::const_ptr poly_list_const_ptr;

	discrete_post_intv_split(double aggregation_type, double clustering) :
		my_aggregation_type(aggregation_type), my_clustering(clustering), my_minbrak_type("gold_desc"),
				my_split_size(0), my_intersection_error(0.0) {
	}
	;
	discrete_post_intv_split(double aggregation_type, double clustering,
				double inters_error) :
			my_aggregation_type(aggregation_type), my_clustering(clustering), my_minbrak_type("gold_desc"),
					my_split_size(0), my_intersection_error(inters_error) {
	}
	;
	discrete_post_intv_split(double aggregation_type, double clustering, std::string minbrak_type,
			double split_size, double inters_error) :
		my_aggregation_type(aggregation_type), my_clustering(clustering), my_minbrak_type(minbrak_type),
				my_split_size(split_size), my_intersection_error(inters_error) {
	}
	;

	virtual ~discrete_post_intv_split() {
	}
	;

	discrete_post_intv_split* clone() {
		return new discrete_post_intv_split<scalar_type> (*this);
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
							"discrete_post_intv_split: continuous_set type "
									+ tname + " not supported as states");
				}
			}
		}
		inst.sinv_ptr = boost::dynamic_pointer_cast<
				const continuous::polyhedron<scalar_type> >(source_inv);
		if (source_inv && !inst.sinv_ptr) {
			throw std::runtime_error(
					"discrete_post_intv_split: continuous_set type not supported as invariant");
		}
		inst.tinv_ptr = boost::dynamic_pointer_cast<
				const continuous::polyhedron<scalar_type> >(target_inv);
		if (target_inv && !inst.tinv_ptr) {
			throw std::runtime_error(
					"discrete_post_intv_split: continuous_set type not supported as invariant");
		}
		inst.guard_poly_ptr = poly_const_ptr();

		if (j.get_guard()) {
			inst.guard_poly_ptr = boost::dynamic_pointer_cast<
					const continuous::polyhedron<scalar_type> >(j.get_guard());
			if (!inst.guard_poly_ptr) {
				throw std::runtime_error(
						"discrete_post_intv_split: guard type not supported");
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
						"discrete_post_sf_inters: transform type not supported\n"
								+ ss.str());
			}
		}
		return inst;
	}
	;
	void post_poly(problem_instance& inst, poly_list_type & polys,
			continuous::constr_polyhedron<scalar_type>& tinv_poly,
			variable_id_set polys_vis) const {
		using namespace continuous;
		using namespace support_function;
		if (inst.reset_trans_ptr) {
			// a pointer for the inputs
			support_function_provider::const_ptr mapped_U;
			math::affine_map<scalar_type> M = *inst.reset_trans_ptr;

			// Remove variables not in the domain of M
			// get the variables from invariant, guard and sfm
			variable_id_set excess_variables = polys_vis;
			set_difference_assign(excess_variables,
					M.domain().get_variable_ids());
			if (!excess_variables.empty()) {
				//LOGGER(DEBUG5,"discrete_post_intv_split::post","quantifying over "+to_string(excess_variables.size())+" variable(s) not in assignment");
				throw std::runtime_error(
						"discrete_post_intv_split:post_poly: excess_variables not empty\n");
				polys.existentially_quantify_variables(excess_variables);
			}

			// If the domain and codomain are not the same,
			// we need to reorder and seperate state from input variables
			if (M.domain() != M.codomain()) {
				// split A into square A1 and A2
				math::affine_map<scalar_type> Mstate;
				math::affine_map<scalar_type> Minput;
				math::separate_states_from_inputs(M, Mstate, Minput);

				M = Mstate;
				//					std::cout << "split dynamics into " << Mstate << " and " << Minput
				//							<< " with input domain " << Minput << std::endl;

				// Define the set of inputs using the invariant.
				// This is the original invariant with both state and nonstate variables
				if (Minput.domain().size() > 0) {
					std::stringstream ss;
					logger::copyfmt_to(ss);
					ss << Minput.domain();
					LOGGER(DEBUG5,"discrete_post_intv_split::post","nondeterministic inputs "+ss.str());
					if (inst.sinv_ptr) {
						ss << inst.sinv_ptr;
						LOGGER(DEBUG5,"discrete_post_intv_split::post","using source invariant "+ss.str());
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
					LOGGER(DEBUG4,"discrete_post_intv_split::post","using outer poly for mapping because of non-regular jump");
				} else if (!singular && mapped_U) {
					LOGGER(DEBUG4,"discrete_post_intv_split::post","using outer poly for mapping because of nondet. inputs");
				} else {
					LOGGER(DEBUG4,"discrete_post_intv_split::post","using outer poly for mapping because of non-regular jump and nondet. inputs");
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
				LOGGERSW(DEBUG4,"discrete_post_intv_split::post","applying affine map to poly in discrete_post_intv_split");
				// transform the polys
				const math::affine_map<scalar_type>& M = *inst.reset_trans_ptr;
				polys = apply_map(polys, M);
			}
		}
		// add constraints with redundancy check because
		// non-redundancy is needed later for template hull
		//polys.add_constraints(*inst.tinv_ptr->get_constraints(), true);
		IFLOGGER(DEBUG5) {
			std::stringstream ss;
			logger::copyfmt_to(ss);
			ss << tinv_poly.get_constraints();
			LOGGER(DEBUG5,"discrete_post_intv_split::post","adding target invariant "+ss.str());
		}
		polys.add_constraints(*tinv_poly.get_constraints(), true);
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
		LOGGERSW(DEBUG1,"discrete_post_intv_split::post","Discrete post with discrete_post_intv_split");

		typedef math::numeric::interval<unsigned int> intv_type;
		typedef std::list<intv_type> intv_list_type;
		typedef constr_polyhedron<scalar_type> poly_type;
		typedef polyhedron_collection<scalar_type> poly_list_type;

		problem_instance inst = get_instance(j, source_inv, target_inv, cset);
		continuous_set_collection ret_collection;

		constr_polyhedron<scalar_type> exit_set_poly, tinv_poly;

		LOGGER(DEBUG2,"discrete_post_intv_split::post","sfm of size "+to_string(inst.sfm_set_ptr->get_size()));
		if (!inst.reset_trans_ptr->is_void()) {

			// if the set is an sfm
			if (inst.sfm_set_ptr) {

				// Initialize the list with all entries of the sfm
				intv_type I; // universe
				I.set_lower(0);
				if (inst.sfm_set_ptr->get_size() > 0)
					I.set_upper(inst.sfm_set_ptr->get_size() - 1);
				else {
					// make the interval finite but empty
					I.set_lower(1);
					I.set_upper(0);
				}

				intv_list_type list;
				if (inst.sfm_set_ptr->get_size() > 0)
					list.push_back(I);

				// Keep track of the variables in the polys
				variable_id_set polys_vis =
						inst.sfm_set_ptr->get_variable_ids();
				set_union_assign(polys_vis, inst.sinv_ptr->get_variable_ids());

				//std::cout << "initially list: " << list << std::endl;

				// @attention We don't intersect with the source invariant because we assume
				// it has already been dealt with.
				// This might change in the future!

				// Find intersecting intervals
				sfm_const_ptr sfm2;

				// need to const_cast because the sfm will be extended
				typename sfm_cont_set<scalar_type>::ptr nonconst_sfm =
						boost::const_pointer_cast<sfm_cont_set<scalar_type> >(
								inst.sfm_set_ptr);

				if (inst.guard_poly_ptr) {
					LOGGER(DEBUG5,"discrete_post_intv_split::post","guard_poly to intersect:"+to_string(inst.guard_poly_ptr));
					list = get_intersection_interval(nonconst_sfm,
							*inst.guard_poly_ptr);
					//std::cout << "getting list from guard: " << list << std::endl;
				} else {
					throw std::runtime_error(
							"discrete_post_intv_split: guard type not supported\n");
				}

				LOGGER(DEBUG2,"discrete_post_intv_split::post","found "+to_string(list.size())+" intervals intersecting with guard");
				/**
				 * DEBUG
				 * print the intersecting intervals
				 */
				IFLOGGER(DEBUG7){
					for (intv_list_type::const_iterator it = list.begin(); it
						!= list.end(); ++it) {
						std::cout << *it << std::endl;
						LOGGER(DEBUG7,"discrete_post_intv_split::post","sfm section intersecting with the guard:"
								+to_string(*it));
					}
				}
				/* Initializations */

				typedef typename std::set<math::vdom_vector<scalar_type>,
						math::numeric::lex_comp_less<scalar_type,
								math::vdom_vector> > vector_set;
				vector_set dir_vector_set;

				/* adding source invariant constraints*/
				constr_polyhedron<scalar_type> guard_poly;
				if (inst.guard_poly_ptr) {
					guard_poly.add_constraints(inst.guard_poly_ptr->get_constraints());
				}
				if (inst.sinv_ptr)
					guard_poly.add_constraints(inst.sinv_ptr->get_constraints());

				typename math::lin_constraint_system<scalar_type>::const_ptr g_cons_sys_ptr;

				/**
				 * Check if the target invariant could be reverse mapped. If so, include the
				 * reverse mapped constraints to the guard constraints list to have
				 * G \cap Source_inv \cap tgt_rev_map_inv
				 */
				std::cout << "Target Invariant:" << inst.tinv_ptr << std::endl;

				LOGGER(DEBUG7,"discrete_post_intv_split::post","guard cons before rev map:" +to_string(guard_poly));
				if(inst.tinv_ptr){
					typename constr_polyhedron<scalar_type>::const_ptr tgt_inv_poly_ptr =
							boost::dynamic_pointer_cast<const constr_polyhedron<scalar_type> >(inst.tinv_ptr);
					if(tgt_inv_poly_ptr){
						math::affine_map<scalar_type> M = *inst.reset_trans_ptr;

						support_function_provider::const_ptr U =
								boost::dynamic_pointer_cast<const support_function_provider>(
										inst.reset_trans_ptr->get_input_set());

						reverse_map_constraints(tgt_inv_poly_ptr, M , U, exit_set_poly, tinv_poly);
						LOGGER(DEBUG5,"discrete_post_intv_split::post","invariant poly reverse mapped to: "+to_string(exit_set_poly));

					}
					else
						throw std::runtime_error("discrete_post: Target location invariant set type not supported\n");
					LOGGER(DEBUG7,"discrete_post_intv_split::post","exit poly:"+to_string(exit_set_poly));
					guard_poly.add_constraints(exit_set_poly.get_constraints());
				}
				LOGGER(DEBUG5,"discrete_post_intv_split::post","guard cons after rev map:" +to_string(guard_poly));
				/**----*/
				guard_poly.remove_redundant_constraints();
				LOGGER(DEBUG5,"discrete_post_intv_split::post","guard constraints after redundancy check : "+to_string(guard_poly));

				if(guard_poly.get_constraints()){
					g_cons_sys_ptr = guard_poly.get_constraints();
				}

				typename math::lin_constraint_system<scalar_type>::ptr
						nonconst_gcons = boost::const_pointer_cast<
								math::lin_constraint_system<scalar_type> >(
								g_cons_sys_ptr);

				continuous_set_ptr c_ptr;

				//collapse guard constraints
				nonconst_gcons->collapse_inequalities();
				LOGGER(DEBUG2,"discrete_post_intv_split::post","guard constraints list after collapsing : "+to_string(nonconst_gcons));
				continuous::constr_polyhedron<scalar_type> cons_poly;
				cons_poly.add_constraints(g_cons_sys_ptr);
				list = get_intersection_interval(nonconst_sfm,cons_poly);
				LOGGER(DEBUG2,"discrete_post_intv_split::post","sfm interval(s) intersecting with the constraint poly(after rev_map): "+to_string(list));

				// construct directions set from the scenario directions choice

				positional_vdomain dom = positional_vdomain(
						nonconst_sfm->get_index_to_variable_id_map());
				typename sfm_cont_set<scalar_type>::vector_list dirs_list =
										 continuous::support_function::direction_chooser::get_directions<scalar_type>(dom);

				for (typename sfm_cont_set<scalar_type>::vector_list::const_iterator
						iter = dirs_list.begin(); iter
						!= dirs_list.end(); ++iter) {
					math::vdom_vector<scalar_type> vdom_vec =
							math::vdom_vector<scalar_type>(dom, *iter);
					dir_vector_set.insert(vdom_vec);
					//std::cout << "direction:" << vdom_vec << std::endl;
				}
				// iterating over the intersecting intervals of the flowpipe.

				for (intv_list_type::const_iterator it = list.begin(); it
						!= list.end(); ++it) {
					std::cout << "treating interval: " << *it << std::endl;
					LOGGER(DEBUG3,"discrete_post_intv_split::post","treating interval of size "+to_string(it->size()));
					// create a new poly collection
					typename poly_list_type::ptr polys_ptr =
							typename poly_list_type::ptr(new poly_list_type());
					poly_list_type& polys = *polys_ptr;

					typename support_function_provider::ptr sf_prov_ptr =
							typename support_function_provider::ptr(
									new continuous::support_function::sfm_section<
											scalar_type>(nonconst_sfm,*it));


					// create a new poly list for containing the outer_polys of
					// only 1 constraint intersection with sfm result.
					typedef std::vector<typename polyhedron<scalar_type>::ptr> poly_vector;
					poly_vector POLYS;

					if (inst.sinv_ptr) {
							if (inst.sinv_ptr->is_empty()) {
							throw std::runtime_error(
								"discrete_post: Source INV set empty\n");
						}
					}

					support_function_provider::ptr g_inters_set_ptr;
					/*
					 * @todo: add guard directions to the outer poly computation
					 */
					bool POLYS_init = false, split_algo;
					typename continuous::support_function::sfm_lin_con_intersection<scalar_type>::ptr sfm_inters_coll;

					// iterating through the guard+src_inv+exit_poly constraints
					for (typename math::lin_constraint_system<
							scalar_type>::const_iterator itr =
							nonconst_gcons->begin(); itr
							!= nonconst_gcons->end(); ++itr) {

						LOGGER(DEBUG2,"discrete_post_intv_split::post","constraint to intersect with sfm section : "+to_string(*itr));
//						std::cout << "POLYS size:" << POLYS.size() << std::endl;

						split_algo = false;
						if (is_intersection_empty(sf_prov_ptr, *itr)) {
							LOGGER(DEBUG2,"discrete_post_intv_split::post","sfm section:"+to_string(*it)+" do not intersect with the constraint\n");
							POLYS = poly_vector();
							break; // Check for next intersecting interval.
						} else if (contains(*itr, sf_prov_ptr)) {
							g_inters_set_ptr = sf_prov_ptr;
							LOGGER(DEBUG2,"discrete_post_intv_split::post","sfm section contained inside the constraint\n");
						} else {
						/*Check for the interval split method algorithm*/
							sfm_inters_coll = typename continuous::support_function::sfm_lin_con_intersection<scalar_type>::ptr(
									new continuous::support_function::sfm_lin_con_intersection<scalar_type>(
											 nonconst_sfm, *it, *itr, my_minbrak_type, my_intersection_error ));
							g_inters_set_ptr = support_function_provider::ptr(sfm_inters_coll);
							split_algo = true;
						}

/*
						std::cout << "Number of directions in the dir_vector_set:" << dir_vector_set.size()
								<< std::endl;
*/

						if(split_algo && sfm_inters_coll){
							// Call the interval split implementation
							poly_vector interm_polys = sfm_inters_coll->get_chull_outer_polys(dir_vector_set,my_split_size);

							if(!POLYS_init){
								POLYS = interm_polys;
								POLYS_init = true;
							}
							else{
								// create a new poly collection
								poly_vector new_polys = poly_vector();
								assert(POLYS.size() == interm_polys.size());

								unsigned int count = 0;
								for(typename poly_vector::const_iterator it1 = POLYS.begin(),
										it2 = interm_polys.begin(); it1 != POLYS.end() && it2!=interm_polys.end(); ++it1, ++it2) {
										typename constr_polyhedron<scalar_type>::ptr cons_poly1 =
												boost::static_pointer_cast<constr_polyhedron<scalar_type> >(*it1);
										typename constr_polyhedron<scalar_type>::ptr cons_poly2 =
												boost::static_pointer_cast<constr_polyhedron<scalar_type> >(*it2);
										cons_poly1->add_constraints(cons_poly2->get_constraints());

										typename polyhedron<scalar_type>::ptr P_ptr = typename polyhedron<scalar_type>::ptr(cons_poly1);
										new_polys.push_back(P_ptr);
								}
								POLYS = new_polys;
							}
						}
						else{
							// create a new poly collection
							poly_vector new_polys = poly_vector();
							poly_vector upd_polys = poly_vector();
							typename support_function_provider::ptr sf_prov_ptr;

							unsigned int l, low = (*it).lower().get_val();
							unsigned int u, up = (*it).upper().get_val();

							assert(my_split_size >= 1);
							size_t intv_size = (*it).size() + 1;
							size_t polys_size = intv_size/my_split_size;
							if(polys_size*my_split_size < intv_size)
								polys_size++;

							for(unsigned int j=0;j<polys_size;j++){
								l = low + j*my_split_size;
								u = l + my_split_size - 1;
								if(u > up)
									u = up;
								scalar_with_infinity<unsigned int> l_scal(l);

								math::numeric::interval<unsigned int> intv(l_scal);
								sf_prov_ptr =
										typename support_function_provider::ptr(
												new continuous::support_function::sfm_section<scalar_type>(nonconst_sfm,intv));
								typename support_function_provider::ptr sf_prov_ptr_1;
								for(unsigned int j = l+1;j<=u;j++){
									intv = math::numeric::interval<unsigned int>(scalar_with_infinity<unsigned int>(j));
									sf_prov_ptr_1 = typename support_function_provider::ptr(
											new continuous::support_function::sfm_section<scalar_type>(nonconst_sfm, intv));

									sf_prov_ptr = typename support_function_provider::ptr(
											new sf_chull<scalar_type>(sf_prov_ptr,sf_prov_ptr_1));
								}

								typename polyhedron<scalar_type>::ptr new_poly =
										typename polyhedron<scalar_type>::ptr(
											new constr_polyhedron<scalar_type> (continuous::compute_outer_poly(*sf_prov_ptr, dir_vector_set)));
								new_polys.push_back(new_poly);
							}

							if(!POLYS_init){
								POLYS = new_polys;
								POLYS_init = true;
							}
							else{

								std::cout << "POLYS size=" << POLYS.size() << std::endl;
								std::cout << "new polys size=" << new_polys.size() << std::endl;
								assert(POLYS.size() == new_polys.size());
								for(typename poly_vector::const_iterator it1 = POLYS.begin(),
									it2 = new_polys.begin(); it1 != POLYS.end() && it2!=new_polys.end(); ++it1, ++it2) {
									typename constr_polyhedron<scalar_type>::ptr cons_poly1 =
											boost::static_pointer_cast<constr_polyhedron<scalar_type> >(*it1);
									typename constr_polyhedron<scalar_type>::ptr cons_poly2 =
											boost::static_pointer_cast<constr_polyhedron<scalar_type> >(*it2);
									cons_poly1->add_constraints(cons_poly2->get_constraints());

									typename polyhedron<scalar_type>::ptr P_ptr = typename polyhedron<scalar_type>::ptr(cons_poly1);
									upd_polys.push_back(P_ptr);
								}
								POLYS = upd_polys;
							}
						}
					}
					/**
					 * Make is poly collection from the list of polys
					 */
					// create a new poly collection
					typename poly_list_type::ptr new_polys_ptr =
							typename poly_list_type::ptr(new poly_list_type());
					poly_list_type& nonred_polys = *new_polys_ptr;

					for(typename poly_vector::const_iterator vec_iter = POLYS.begin(); vec_iter != POLYS.end(); vec_iter++){
						nonred_polys.insert(*vec_iter);
					}
					//-------
					if (!nonred_polys.is_empty()) {
						set_union_assign(polys_vis,
								nonred_polys.get_variable_ids());
						polys.insert(nonred_polys);
					}
/*
					std::cout << "discrete_post_intv_split: polys after intersection with guard constraint follows, size=" << polys.size();
					for (typename polyhedron_collection<scalar_type>::const_iterator
													it = polys.begin(); it != polys.end(); ++it) {
						std::cout << *it << std::endl;

					}
*/
					polys.remove_empty();
					/* DEBUG */
					if (polys.is_empty()) {
						LOGGER(DEBUG3,"discrete_post_intv_split:post:","polys list empty" );
					}

					IFLOGGER(DEBUG7) {
						for (typename polyhedron_collection<scalar_type>::const_iterator
								it = polys.begin(); it != polys.end(); ++it) {

							if ((*it)->is_universe()) {
								LOGGER(DEBUG7,"discrete_post_intv_split:post:","Intersection Poly is Universe!");
							}
							//ret_collection.insert(*it);
						}
					}

					post_poly(inst, polys, tinv_poly, polys_vis);

					// either pass on the convex hull of the polys
					// or the template hull...

					if (polys.size() > 1 && my_clustering > 0.0
							&& my_aggregation_type != 1) {
						// there's no need to do this if aggregation is thull without size limit
						LOGGERSW(DEBUG3,"discrete_post_intv_split::post","clustering");
						LOGGER_ATTACH(DEBUG3,"discrete_post_intv_split::post"," "+to_string(polys.size()) + " polyhedra using template hull, ",logger::get_last_id());
						double err = my_clustering;
						polys
								= constraint_hull_generator<scalar_type> (
										*polys_ptr).compute_constraint_hull_up_to_error(
										err);

						LOGGER_ATTACH(DEBUG3,"discrete_post_intv_split::post"," obtained "+to_string(polys.size()) + " polyhedral clusters, ",logger::get_last_id());
					}
					if (!polys.is_empty()) {
						if (my_aggregation_type == 0 || double(polys.size())
								<= my_aggregation_type) {
							// use convex hull, i.e., insert the poly_collection to be used as sf_provider

							// if it's just one set, unwrap
							if (polys.size() == 1) {
								LOGGER(DEBUG3,"discrete_post_intv_split::post","adding " + to_string(polys.size()) + " polyhedron to result");
								ret_collection.insert(*polys.begin());
							} else {
								LOGGER(DEBUG3,"discrete_post_intv_split::post","adding convex hull of " + to_string(polys.size()) + " polyhedra to result");
								ret_collection.insert(polys_ptr);
							}
						} else if (my_aggregation_type > 0) {
							// use template hull
							LOGGERSW(DEBUG3,"discrete_post_intv_split::post","adding polyhedral template hull");
							LOGGER_ATTACH(DEBUG3,"discrete_post_intv_split::post"," of "+to_string(polys.size()) + " polyhedra to result, ",logger::get_last_id());
							ret_collection.insert(
									constraint_hull_generator<scalar_type> (
											*polys_ptr).compute_constraint_hull());
						} else {
							// treat each convex set separately
							//  -- creating a new symbolic state for each one
							LOGGER(DEBUG3,"discrete_post_intv_split::post","adding " + to_string(polys.size()) + " separate polyhedra to result");
							for (typename polyhedron_collection<scalar_type>::const_iterator
									it = polys.begin(); it != polys.end(); ++it) {
								ret_collection.insert(*it);
							}
						}
					} else {
						LOGGER(DEBUG3,"discrete_post_intv_split::post","result of jump is empty");
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

				LOGGER(DEBUG2,"discrete_post_intv_split::post","computing discrete post of polyhedron");

				// create a new poly collection
				// and add the poly to it
				poly_list_type polys;
				if (inst.cset_poly_coll_ptr) {
					polys.insert(*inst.cset_poly_coll_ptr);
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

				post_poly(inst, polys, tinv_poly, polys_vis);

				for (typename polyhedron_collection<scalar_type>::const_iterator
						it = polys.begin(); it != polys.end(); ++it) {
					ret_collection.insert(*it);
				}
			}
		}
		return ret_collection;
}
;

private:
	double my_aggregation_type;
	double my_clustering;
	double my_intersection_error;
	size_t my_split_size;
	std::string my_minbrak_type;
};

}

#endif

