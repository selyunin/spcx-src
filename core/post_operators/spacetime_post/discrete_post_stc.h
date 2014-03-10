/*
 * discrete_post_stc.h
 *
 *  Created on: Nov 3, 2012
 *      Author: notroot
 */

#ifndef DISCRETE_POST_STC_H_
#define DISCRETE_POST_STC_H_

#include <typeinfo>

#include "utility/logger_stopwatch.h"

#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "core/hybrid_automata/transition.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/post_operators/discrete_post.h"
#include "core/continuous/support_function/spacetime_flowpipe.h"
#include "core/continuous/support_function/sf_derived/sf_set_operations.h"
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
class discrete_post_stc: public discrete_post {
public:
	typedef boost::shared_ptr<discrete_post_stc<scalar_type> > ptr;
	typedef boost::shared_ptr<const discrete_post_stc<scalar_type> > const_ptr;

	typedef continuous::spacetime_flowpipe<scalar_type> flowpipe;
	typedef typename flowpipe::ptr flowpipe_ptr;
	typedef typename flowpipe::const_ptr flowpipe_const_ptr;
	typedef typename flowpipe::cut_point_method cut_point_method;

	typedef typename flowpipe::time_interval time_interval;
	typedef std::vector<time_interval> intv_list_type;

	typedef typename continuous::constr_polyhedron<scalar_type> cpoly_type;
	typedef typename continuous::polyhedron<scalar_type> poly_type;
	typedef typename continuous::polyhedron<scalar_type>::const_ptr poly_const_ptr;
	typedef typename continuous::polyhedron<scalar_type>::ptr poly_ptr;
	typedef continuous::polyhedron_collection<scalar_type> poly_list_type;
	typedef typename poly_list_type::const_ptr poly_list_const_ptr;

	typedef math::affine_map<scalar_type> affine_map;
	typedef boost::shared_ptr<const continuous::support_function_provider> supp_const_ptr;
	typedef continuous::reset_affine_transform<scalar_type> reset_type;
	typedef boost::shared_ptr<const reset_type > reset_type_const_ptr;

	typedef typename flowpipe::error_type error_type;

	discrete_post_stc(const cut_point_method& m) {
		my_cut_point_method = m;
	}
	;

	virtual ~discrete_post_stc() {
	}
	;

	discrete_post_stc* clone() {
		return new discrete_post_stc<scalar_type>(*this);
	}
	;

	/** Obtains a polyhedron (in space-time) from a spacetime_flowpipe
	 *
	 * The quantification over the time variable and the last discrete
	 * jump are taken care of by mapping. */
	static void get_poly_from_flowpipe(flowpipe_const_ptr& cflow,
			error_type err) {
		using namespace continuous;
		using namespace support_function;

		throw std::runtime_error("missing implementation");

		// -------------------------------------------------------------
		// if cset is a flowpipe piece, compute its outer polytope and map it
		// -------------------------------------------------------------
		// get a single outer polytope
		typename flowpipe::cut_point_method m;
		m.type = flowpipe::cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES;
		m.piece_count = 1;
		m.approx_error = err;
		polyhedron_collection<scalar_type> polys =
				cflow->compute_outer_polyhedra(m);
		// lock outer approximation so that containment checking can be performed
		flowpipe* nonconst_cflow = const_cast<flowpipe*>(cflow.get());
		nonconst_cflow->set_refinement_locked();

		assert(polys.size() == 1);

		// Now map the polyhedra (in space-time)
		if (!math::definitely((*polys.begin())->is_empty())) {
//std::cout << std::endl << "Poly: " << **polys.begin() << std::endl;
			support_function_provider::const_ptr mapped_set;
			if (cflow->get_map()) {
//std::cout << "Map: " << *cflow->get_map() << std::endl;
				mapped_set = map_and_bloat(*polys.begin(), *cflow->get_map(),
						cflow->get_bloating());
			} else {
				mapped_set = bloat<scalar_type>(*polys.begin(),
						cflow->get_bloating());
			}
			IFLOGGER(DEBUG5) {
				hyperbox<scalar_type> bb = compute_bounding_box<scalar_type>(
						*mapped_set);
				LOGGER_OS(DEBUG5, __FUNCTION__) << "bounding box: " << bb
						<< std::endl;
			}
			cflow = mapped_set;
		} else {
			cflow = flowpipe_const_ptr();
		}
	}
	;

	/** Returns a sequence of flowpipe segments, each of which intersects with the guard
	 *
	 * The segmentation is done such that the number of segments is minimized within
	 * the error bounds.
	 *
	 * The convexification is done here as well, since the segmentation is tightened
	 * after the convexification.
	 * */
	std::vector<flowpipe_ptr> intersect_flowpipe(const flowpipe& flow,
			const math::lin_constraint_system<double>& cons) const {
		using namespace continuous;

		std::vector<flowpipe_ptr> flowpipe_segments;

		// Initialize the list with the entire domain of the flowpipe
		intv_list_type list(1, time_interval::whole());
		intv_list_type narrow_list(1, time_interval::whole());

		// Find intersecting intervals
		LOGGERSWOC(DEBUG5, __FUNCTION__, "Getting guard interval");
		// need to const_cast because evolutions will be added
		flowpipe& nonconst_flow = const_cast<flowpipe&>(flow);

		/** Obtain minimal number of segments by first intersecting with
		 *  bloated lower bound, then with actual upper bound and
		 *  combining ithe intervals.
		 */
		bool used_narrow_list = true;
		narrow_list = nonconst_flow.get_maybe_satisfying_subdomains(cons,
				my_cut_point_method.approx_error, false);

		if (narrow_list.size() > 1) {
			list = nonconst_flow.get_maybe_satisfying_subdomains(cons,
					my_cut_point_method.approx_error, true);
			// delimit the list by min and max of narrow_list
			list = plif::tighten_subintervals(list, narrow_list);
			LOGGER(DEBUG3, "discrete_post_stc::post",
					"used convex hull to reduce from "
							+ to_string(narrow_list.size()) + " to "
							+ to_string(list.size())
							+ " intersecting time intervals");
			used_narrow_list = false;
		} else {
			list = narrow_list;
			used_narrow_list = true;
		}

		LOGGER(DEBUG2, __FUNCTION__,
				"found " + to_string(list.size())
						+ " intervals intersecting with guard");

		/**
		 * Creating a segment for each interval intersecting with the guard
		 *  */
		for (typename intv_list_type::const_iterator it = list.begin();
				it != list.end(); ++it) {
			LOGGER(DEBUG3, __FUNCTION__,
					"creating segment for time interval [" + to_string(it->lower()) + ","
							+ to_string(it->upper()) + "]");

			// make a copy of the flowpipe and restrict it to the current interval
			flowpipe_ptr guard_flow(flow.clone());
			{
				LOGGERSWOC(DEBUG5, __FUNCTION__, "Intersection with guard");

				guard_flow->restrict_to_subdomain(*it);
				// intersect with the guard
				if (!cons.empty()) {
					// restrict_evolutions = true
					guard_flow->add_outer_constraints(cons, true);
				}
			}

			/** Convexify */
			std::vector<spacetime_flowpipe<scalar_type> > convex_guard_pieces;
			convex_guard_pieces = guard_flow->convexify(my_cut_point_method);
			if (convex_guard_pieces.size() > 1) {
				LOGGER(DEBUG4, __FUNCTION__,
						"convexification resulted in "
								+ to_string(convex_guard_pieces.size())
								+ " pieces");
			}

			// add to result
			for (typename std::vector<spacetime_flowpipe<scalar_type> >::const_iterator it =
					convex_guard_pieces.begin();
					it != convex_guard_pieces.end(); ++it) {
				typename spacetime_flowpipe<scalar_type>::ptr p(
						new spacetime_flowpipe<scalar_type>(*it));

				/** Restrict again after convexification */
				if (!used_narrow_list) {
					intv_list_type tightened_list(1, time_interval::whole());
					tightened_list = p->get_maybe_satisfying_subdomains(cons,
							my_cut_point_method.approx_error, false);
					time_interval tightened_intv = plif::convex_hull(tightened_list);
					if (tightened_intv.is_empty())
						throw std::runtime_error("unexpected empty interval");
					p->restrict_to_subdomain(tightened_intv);
				}

				// lock outer approximation so that containment checking can be performed
				// -> do that at the latest possible occasion, i.e., when the outer poly is actually used
				//p->set_refinement_locked();

//				std::cout << "got flowpipe segment: " << *p;

				flowpipe_segments.push_back(p);
			}
		}

		return flowpipe_segments;
	}

	struct problem_instance {
		flowpipe_const_ptr sflow_ptr;
		poly_list_const_ptr cset_poly_coll_ptr;
		poly_const_ptr cset_poly_ptr;
		poly_const_ptr sinv_ptr;
		poly_const_ptr tinv_ptr;
		poly_const_ptr guard_poly_ptr;
		affine_map M;
		supp_const_ptr W;
	};

	problem_instance get_instance(const jump_constraints& j,
			continuous::continuous_set::const_ptr source_inv,
			continuous::continuous_set::const_ptr target_inv,
			continuous::continuous_set::const_ptr cset) const {
		using namespace continuous;
		problem_instance inst;

		inst.guard_poly_ptr = poly_const_ptr();
		if (j.get_guard()) {
			inst.guard_poly_ptr = boost::dynamic_pointer_cast<
					const polyhedron<scalar_type> >(j.get_guard());
			if (!inst.guard_poly_ptr) {
				throw std::runtime_error(
						"discrete_post_stc: guard type not supported");
			}
		}
		inst.M = affine_map();
		if (j.get_transform()) {
			reset_type_const_ptr trans_p =
					boost::dynamic_pointer_cast<const reset_type>(j.get_transform());
			if (!trans_p) {
				std::stringstream ss;
				logger::copyfmt_to(ss);
				ss << "offending transform:" << j.get_transform();
				throw std::runtime_error(
						"discrete_post_stc: transform type not supported\n"
								+ ss.str());
			} else {
				inst.M = *trans_p;
				inst.W = boost::static_pointer_cast<
									const support_function_provider>(
									trans_p->get_input_set());
			}
		}

		inst.sflow_ptr = boost::dynamic_pointer_cast<const flowpipe>(cset);
		// if it's a flowpipe check first whether it's to be
		// mapped and bloated
		if (inst.sflow_ptr) {
			// If to be bloated, we need to use the outer poly
			if (inst.sflow_ptr->get_bloating()) {
				// The flowpipe is mapped with Y = MX+b+V

				// map is nondeterministic, so use outer poly
//				get_poly_from_flowpipe(inst.sflow_ptr,my_cut_point_method.approx_error);

				// Note: this case should be taken care of by continuous_post (for simplicity)
				throw std::runtime_error("missing implementation");
			} else if (inst.sflow_ptr->get_map()) {
				LOGGER(DEBUG5, __FUNCTION__,
						"joining transitions after skipping time elapse");

				// The flowpipe is mapped with Y = KX+b, but
				// there is no bloating (V=0)
				const affine_map& K = *inst.sflow_ptr->get_map();

				// the map is deterministic, so concatenate
				// 1. pull up the guard intersection to before M
				if (inst.guard_poly_ptr) {
					LOGGER_OS(DEBUG7,__FUNCTION__) << "pre-image of guard: " << inst.guard_poly_ptr << " with map " << K;
					const constr_polyhedron<scalar_type>& G =
							static_cast<const constr_polyhedron<scalar_type>&>(*inst.guard_poly_ptr);
					inst.guard_poly_ptr = poly_ptr(
							new constr_polyhedron<scalar_type>(
									reverse_map(G, K)));
					LOGGER_OS(DEBUG7,__FUNCTION__) << "pre-mapped guard: " << inst.guard_poly_ptr;
				}
				// 2. Concatenate maps
				if (inst.M.is_empty()) {
					// new map : K*X + b + W
					inst.M = K;
				} else {
					// new map : M*K*X + M*v + (w+W)
					inst.M = concatenate(inst.M, K);
					// note: W stays
				}
				LOGGER_OS(DEBUG7,__FUNCTION__)
						<< "got map with domain: " << inst.M.domain() << ", "
						<< "affine map: " << inst.M << std::flush;

				/** @note We should really set the map of inst.sflow_ptr to null
				 * because it now has been fused with the new transition's map. */
				typename spacetime_flowpipe<scalar_type>::ptr p(
						new spacetime_flowpipe<scalar_type>(*inst.sflow_ptr));
				p->set_map(affine_map());
				inst.sflow_ptr = p;
				// std::cout << "got flowpipe: " << *inst.sflow_ptr;
			}
		}

		if (!inst.sflow_ptr) {
			inst.cset_poly_coll_ptr = boost::dynamic_pointer_cast<
					const poly_list_type>(cset);
			if (!inst.cset_poly_coll_ptr) {
				inst.cset_poly_ptr = boost::dynamic_pointer_cast<
						const polyhedron<scalar_type> >(cset);
				if (!inst.cset_poly_ptr) {
					std::string tname = typeid(*cset.get()).name();
					throw std::runtime_error(
							"discrete_post_stc: continuous_set type " + tname
									+ " not supported as states");
				}
			}
		}

		inst.sinv_ptr = boost::dynamic_pointer_cast<
				const polyhedron<scalar_type> >(source_inv);
		if (source_inv && !inst.sinv_ptr) {
			throw std::runtime_error(
					"discrete_post_stc: continuous_set type not supported as invariant");
		}
		inst.tinv_ptr = boost::dynamic_pointer_cast<
				const polyhedron<scalar_type> >(target_inv);
		if (target_inv && !inst.tinv_ptr) {
			throw std::runtime_error(
					"discrete_post_stc: continuous_set type not supported as invariant");
		}
		return inst;
	}
	;

	void post_poly(problem_instance& inst, poly_list_type & polys,
			variable_id_set polys_vis) const {

		LOGGERSWOC(DEBUG4, __FUNCTION__,
				"Computing discrete successors of polyhedra");

		using namespace continuous;
		using namespace support_function;
		if (!inst.M.is_empty()) {
			// a pointer for the inputs
			support_function_provider::const_ptr mapped_U;
			affine_map& M = inst.M;

			// Remove variables not in the domain of M
			// get the variables from invariant, guard and sfm
			variable_id_set excess_variables = polys_vis;
			set_difference_assign(excess_variables,
					M.domain().get_variable_ids());
			if (!excess_variables.empty()) {
				if (quantify_unused_jump_variables) {
					// Use existential quantification to eliminate variables
					LOGGER(DEBUG5, "discrete_post_stc::post",
							"quantifying over "+to_string(excess_variables.size())+" variable(s) not in assignment");

					// additional log output
					IFLOGGER(DEBUG7) {
						std::stringstream ss;
						logger::copyfmt_to(ss);
						print_variable_id_set(ss, polys_vis);
						LOGGER_OS(DEBUG7,"discrete_post_stc::post")
								<< "poly variables: " << ss.str() << ", ";
						ss.str("");
						print_variable_id_set(ss,
								M.domain().get_variable_ids());
						LOGGER_OS(DEBUG7,"discrete_post_stc::post")
								<< "map domain: " << ss.str() << ", "
								<< "affine map: " << M << std::flush;
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

					positional_vdomain new_dom = M.domain();
					new_dom.add_variables(
							create_variable_set(excess_variables));
					M.reorder(M.codomain(), new_dom);

					IFLOGGER(DEBUG7) {
						// additional log output
						LOGGER_OS(DEBUG7,"discrete_post_stc::post")
								<< "Changed domain of jump map to include variables of continuous set: ";
						std::stringstream ss;
						logger::copyfmt_to(ss);
						print_variable_id_set(ss, polys_vis);
						LOGGER_OS(DEBUG7,"discrete_post_stc::post") << ss.str()
								<< ". New matrix:" << M << std::endl
								<< "new codomain:" << M.codomain()
								<< ", new domain:" << M.domain() << std::endl;
					}
				}
			}

			/** Transform the map and the source invariant into state map and input set */
			// compute_map_and_inputs(M, inst.sinv_ptr, mapped_U);
			mapped_U = inst.W;


			// If A is singular or there are nondet. inputs,
			// map the support function and the obtain an outer poly approx.
			// Otherwise, map the polys directly (exact solution).
			scalar_type det = math::matrix_determinant(M.get_A().get_matrix());
			bool singular = math::numeric::is_MEQ(det / scalar_type(1000),
					scalar_type(0));
			if (singular || mapped_U) {
				if (singular && !mapped_U) {
					LOGGER(DEBUG4, "discrete_post_stc::post",
							"using outer poly for mapping because of non-regular jump");
				} else if (!singular && mapped_U) {
					LOGGER(DEBUG4, "discrete_post_stc::post",
							"using outer poly for mapping because of nondet. inputs");
				} else {
					LOGGER(DEBUG4, "discrete_post_stc::post",
							"using outer poly for mapping because of non-regular jump and nondet. inputs");
				}
				// ... use sf to transform and get outer poly

				// get the directions for the outer polys
				std::list<math::vector<scalar_type> > directions;
				choose_directions(M.codomain(), directions); // R is the number of directions.
				typename sf_set<scalar_type>::vector_set sf_directions;
				for (typename std::list<math::vector<scalar_type> >::const_iterator it =
						directions.begin(); it != directions.end(); ++it) {
					typename sf_set<scalar_type>::vector_type sf_vec(
							M.codomain(), *it);
					sf_directions.insert(sf_vec);
				}

				// make a new list of polys with the result
				poly_list_type mapped_polys;
				for (typename poly_list_type::const_iterator it = polys.begin();
						it != polys.end(); ++it) {
					// the result is the minkowski sum of the mapped states and the mapped inputs
					typename sf_set<scalar_type>::ptr mapped_sf = map_and_bloat(
							*it, M, mapped_U);

					// get the outer poly of the mapped sf object
					typename sf_set<scalar_type>::poly_ptr mapped_poly =
							mapped_sf->outer_poly(sf_directions);
					// add the poly to the result
					mapped_polys.insert(mapped_poly);
				}

				// replace polys with the mapped polys
				polys.swap(mapped_polys);
			} else {
				LOGGERSW(DEBUG4, "discrete_post_stc::post",
						"applying affine map to poly in discrete_post_stc");
				// transform the polys
				polys = apply_map(polys, inst.M);
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
				LOGGER(DEBUG5, "discrete_post_stc::post",
						"adding target invariant "+ss.str());
			}
			polys.add_constraints(*inst.tinv_ptr->get_constraints(), true);
			polys.remove_empty();
		}
		polys.remove_empty();
	}

	/** Compute the post image of a spacetime_flowpipe by taking its outer polyhedral approximation
	 */
	continuous::continuous_set_collection post(const jump_constraints& j,
			continuous::continuous_set::const_ptr source_inv,
			continuous::continuous_set::const_ptr target_inv,
			continuous::continuous_set::const_ptr cset) const {
		using namespace continuous;
		using namespace support_function;
		LOGGERSW(DEBUG1, "discrete_post_stc::post", "Discrete post");

		problem_instance inst = get_instance(j, source_inv, target_inv, cset);
		continuous_set_collection ret_collection;

		if (!inst.M.is_void()) {

			// If the set is an spacetime_flowpipe
			if (inst.sflow_ptr && (!inst.guard_poly_ptr || !math::definitely(inst.guard_poly_ptr->is_empty()))) {
				/**
				 * Restrict flowpipe to intervals intersecting with the guard
				 *  */

				LOGGERSW(DEBUG2, __FUNCTION__,
						"computing discrete post of spacetime_flowpipe");
//				LOGGER_ATTACH(DEBUG2, "discrete_post_stc::post",
//						" of size "+to_string(inst.sflow_ptr->get_size()),
//						logger::get_last_id());

				math::lin_constraint_system<double> guard_cons;
				if (inst.guard_poly_ptr) {
					guard_cons = *inst.guard_poly_ptr->get_constraints();
				}

				/** Compute the segments that
				 * a) intersect with the guard
				 * b) are convex.
				 */
				std::vector<flowpipe_ptr> flowpipe_segments;
				flowpipe_segments = intersect_flowpipe(*inst.sflow_ptr,
						guard_cons);

				/** Transform the map and the source invariant into state map and input set */
				affine_map M = inst.M;

				/** Extend the map so that it quantifies over spacetime */
				positional_vdomain dom = M.domain();
				positional_vdomain domt = inst.sflow_ptr->domain();
				domt.add_variable(inst.sflow_ptr->get_time_variable());

//					std::cout << std::endl << "Map Domain:" << M.domain() << std::endl;
//					std::cout << "Map: " << M << std::endl;

				try {
					if (true) {
						affine_map M_proj;
						M_proj = affine_map::projection_map(M.domain(), domt);
						M = concatenate(M, M_proj);
					} else {
						M.reorder(M.codomain(), domt);
					}
				} catch (basic_exception& e) {
					LOGGER_OS(LOW,__FUNCTION__)
							<< "trying to reorder domain of map from " << dom
							<< " to " << domt;
					LOGGER_OS(LOW,__FUNCTION__) << "Map = " << M;
					LOGGER_OS(LOW,__FUNCTION__) << "Flowpipe = " << *inst.sflow_ptr;
					throw basic_exception(
							"Could not find one of the jump variables in the flowpipe description. Please declare algebraic variables as uncontrolled.",
							e);
				}
//std::cout << std::endl << "Proj: " << M_proj << " Map: " << M << std::endl;

				/** Map each segment and add it to the result */
				for (typename std::vector<flowpipe_ptr>::iterator it =
						flowpipe_segments.begin();
						it != flowpipe_segments.end(); ++it) {
					flowpipe_ptr& p = *it;
					// assign map and bloat set to guard flow
					p->set_map(M);
					p->set_bloating(inst.W);

					ret_collection.insert(p);
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

				LOGGERSW(DEBUG2, "discrete_post_stc::post",
						"computing discrete post of polyhedron");

				// create a new poly collection
				// and add the poly to it
				poly_list_type polys;
				if (inst.cset_poly_coll_ptr) {
					// clone the polys
					for (typename polyhedron_collection<scalar_type>::const_iterator it =
							inst.cset_poly_coll_ptr->begin();
							it != inst.cset_poly_coll_ptr->end(); ++it) {
						polys.insert(
								typename polyhedron<scalar_type>::ptr(
										(*it)->clone()));
					}
				} else {
					polys.insert(
							typename polyhedron<scalar_type>::ptr(
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

				for (typename polyhedron_collection<scalar_type>::const_iterator it =
						polys.begin(); it != polys.end(); ++it) {
					// map each poly
					ret_collection.insert(*it);
				}
			}
		} else {
			// null transform pointer
			throw basic_exception(
					"internal error: unexpected null pointer for discrete transform");
		}

		//		else {
		//				// the transform is empty
		//				continuous::continuous_set_ptr ret_set
		//				= continuous::continuous_set_ptr(
		//						new continuous::constr_polyhedron<scalar_type>(
		//								continuous::constr_polyhedron<scalar_type>::empty_poly()));
		//				ret_collection.insert(ret_set);
		//			}

		if (ret_collection.size() > 1) {
			LOGGER(HIGH, __FUNCTION__,
					"discrete transition produced "
							+ to_string(ret_collection.size()) + " pieces");
		}
		return ret_collection;
	}
	;

	/** Overloard of the standard transition post that add the states also to passed list */
	virtual void add_post_states_trans(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const transition_id& trans, const symbolic_state_ptr& sstate) const {
		// call the base operator
		discrete_post::add_post_states_trans(aut,passed_result_set,waiting_result_set,trans,sstate);
		// copy wait states into passed
		passed_result_set->union_assign(waiting_result_set);
	}

private:
	cut_point_method my_cut_point_method;

public:
	static bool quantify_unused_jump_variables;
};

template<class scalar_type>
bool discrete_post_stc<scalar_type>::quantify_unused_jump_variables = true;

}

#endif /* DISCRETE_POST_STC_H_ */
