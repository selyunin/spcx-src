/*
 * continuous_post_stc.h
 *
 *  Created on: Nov 3, 2012
 *      Author: notroot
 */

#ifndef CONTINUOUS_POST_STC_H_
#define CONTINUOUS_POST_STC_H_

#include "utility/logger_stopwatch.h"

#include "core/post_operators/continuous_post.h"

#include <typeinfo>
#include <queue>

//#include "../abstract_framework/global_types.h"
//#include "../utility/shared_ptr_output.h"
//#include "../valuation_function/node_print_visitor.h"
#include "utility/basic_warning.h"

#include "core/continuous/continuous_dynamics/ode_affine_dynamics.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/continuous/continuous_set_collection.h"
#include "core/hybrid_automata/location.h"
#include "core/continuous/polyhedra/polyhedron.h"
#include "core/continuous/polyhedra/hyperbox/bounding_box_utility.h"
#include "core/continuous/support_function/spacetime_flowpipe.h"
#include "core/continuous/support_function/template_directions/choose_directions.h"
#include "core/continuous/support_function_provider.h"
#include "core/continuous/support_function_provider_utility.h"
#include "core/continuous/support_function/sf_derived/sf_set_operations.h"

#include "core/post_operators/dynamics_preprocessing.h"

// interface to NLtoolbox
#include "math/numeric/boost_interval_utility.h"
#include "core/post_operators/nonlinear_post/affine_interpolation.h"
// for computing bounding box
#include "core/continuous/support_function/spacetime_flowpipe_utility.h"

// for debug output
#include "io/GEN_format/GEN_formatter.h"

namespace hybrid_automata {

/** \brief A continuous post operator for ode_affine_dynamics based
 * on spacetime_plif
 *
 * It is applicable under the following conditions:
 * - the dynamics are of type ode_affine_dynamics,
 * - the continuous_set class is a support_function_provider.
 */

template<class scalar_type>
class continuous_post_stc: public continuous_post {
public:
	typedef boost::shared_ptr<continuous_post_stc> ptr;
	typedef boost::shared_ptr<const continuous_post_stc> const_ptr;

	typedef math::affine_map<scalar_type> affine_map;
	typedef math::vdom_matrix<scalar_type> matrix_type;

	typedef typename continuous::spacetime_flowpipe<scalar_type> flowpipe;
	typedef typename flowpipe::error_type error_type;
	typedef typename continuous::polyhedron<scalar_type>::ptr poly_ptr;

	struct time_elapse_parameters {
		time_elapse_parameters() :
				linearization_error_factor(10), linearization_upscale_factor(2), domain_distance_factor(
						10), refine_linearization_error(true), refine_linearization_with_reach(false), default_dwell_time_factor(
						0.1), nonlinear_flowcut(false), split_linearization(true) {
		}
		;
		/** Error used for affine reachability computation*/
		typename continuous::spacetime_flowpipe<scalar_type>::error_type affine_error;
		/** The linearization error is obtained by multiplying the affine_error with this factor
		 * @attention Currently, only the absolute error is used. */
		scalar_type linearization_error_factor;
		/** If the linearization fails, augment error bound by this factor */
		scalar_type linearization_upscale_factor;
		/** If the initial set is closer than this distance*affine_error.abs() to the domain, increase linearization error. */
		scalar_type domain_distance_factor;
		/** Refine linearization error with a-posteriori measurements */
		bool refine_linearization_error;
		/** Refine linearization error with reachable states and repeat reach */
		bool refine_linearization_with_reach;
		/** Default dwell time limit
		 *
		 * If the dwell time estimation doesn't give a result, use time_horizon*factor
		 */
		scalar_type default_dwell_time_factor;
		/** Use domain intersection in flow rather than last set in domain */
		bool nonlinear_flowcut;
		/** Split set if linearization error too large. If false, the error bound is increased instead. */
		bool split_linearization;
	};

	continuous_post_stc(const time_elapse_parameters& param) :
			my_param(param) {
	}
	;

	virtual ~continuous_post_stc() {
	}
	;

	virtual continuous_post_stc* clone() {
		return new continuous_post_stc(*this);
	}
	;

	/** Intersect Flowpipe with Invariant */
	virtual void restrict_to_invariant(typename flowpipe::ptr& sfp_ptr,
			const math::lin_constraint_system<scalar_type>& inv_constraints, typename flowpipe::error_type err) const {
		using namespace math;
		using namespace continuous;
		using namespace support_function;

		LOGGER(DEBUG3, "continuous_post_stc::restrict_to_invariant",
				"intersecting flowpipe with invariant");

		// Intersect the result with the invariant
		// The result is not empty because the initial set is not empty.
		if (!inv_constraints.empty()) {
			/** We now assume the constraints are already reduced to state variables */
//			const positional_vdomain& dom =
//					sfp_ptr->get_flowpipe_description().dyn.codomain();
//			lin_constraint_system<double> cons =
//					compute_state_variable_constraints<scalar_type>(dom, inv_constraints);

			// restrict domain to when invariant is satisfied
			// @todo find domain end by iteratively refining the error
			sfp_ptr->restrict_to_maybe_sat_prefix(inv_constraints, err);
			LOGGER_OS(DEBUG7,__FUNCTION__) << "restricted to sat prefix of " << inv_constraints << ", resulting in time domain " << sfp_ptr->get_time_domain() << std::endl;
			sfp_ptr->add_outer_constraints(inv_constraints);
		}
	}

	/** Restrict Flowpipe to time domain in which it satisfies constraints */
	virtual void restrict_to_sat_domain(typename flowpipe::ptr& sfp_ptr,
			const math::lin_constraint_system<scalar_type>& sat_constraints, typename flowpipe::error_type err) const {
		using namespace math;
		using namespace continuous;
		using namespace support_function;

		LOGGER(DEBUG3, "continuous_post_stc::restrict_to_invariant",
				"intersecting flowpipe with invariant");

		// Intersect the result with the invariant
		// The result is not empty because the initial set is not empty.
		if (!sat_constraints.empty()) {
			const positional_vdomain& dom =
					sfp_ptr->get_flowpipe_description().dyn.codomain();
			// We assume that we have only state constraints, so the following projection is not necessary
			// lin_constraint_system<double> cons = compute_state_variable_constraints(dom,inv_constraints);

			// restrict domain to when invariant is satisfied
			sfp_ptr->restrict_to_definitely_sat_prefix(sat_constraints, err);
		}
	}

	/** Add User Template Directions */
	virtual void add_user_template_directions(
			typename flowpipe::ptr& sfp_ptr) const {
		using namespace math;
		using namespace continuous;
		using namespace support_function;

		// error
		typename flowpipe::error_type err = my_param.affine_error;

		// add template directions
		const positional_vdomain& dom =
				sfp_ptr->get_flowpipe_description().dyn.codomain();

		LOGGER(DEBUG3, "continuous_post_stc::add_user_template_directions",
				"adding user template directions for "+to_string(dom.size())+" state variables");

		typedef std::list<math::vector<double> > direction_list_type;
		direction_list_type dir_list =
				direction_chooser::get_directions<double>(dom);
		typename flowpipe::direction d(dom);
		for (direction_list_type::const_iterator it = dir_list.begin();
				it != dir_list.end(); ++it) {
			d.set_vector(*it);
			sfp_ptr->get_or_compute_evolution(d, err);
		}
	}

	/** Creates a flowpipe object for given initial states, invariant, affine dynamics, starting at time_offset */
	virtual typename flowpipe::ptr create_flowpipe_object(
			const continuous::support_function_provider::const_ptr& csup,
			const math::lin_constraint_system<scalar_type>& inv_constraints,
			const affine_map& M,
			const continuous::support_function_provider::const_ptr& U,
			double time_offset, double time_horizon) const {
		using namespace math;
		using namespace continuous;
		using namespace support_function;

		// describe the problem
		typename flowpipe::flowpipe_description aff_descr;

		aff_descr.X0 = csup;
		aff_descr.dyn = continuous::ode_affine_dynamics<scalar_type>(M,U);
		// add invariant to boundary conditions (currently not used)
		if (!inv_constraints.empty()) {
			// pass a copy of the invariant constraints
			typename math::lin_constraint_system<scalar_type>::ptr cons_ptr(
					new math::lin_constraint_system<scalar_type>(
							inv_constraints));
			aff_descr.V = cons_ptr;
		}

		const positional_vdomain& dom = aff_descr.dyn.codomain();

		// compute the time interval over which the flowpipe ranges
		typename flowpipe::time_interval tdom(time_offset, time_horizon);
		typename flowpipe::ptr sfp_ptr(new flowpipe(aff_descr, tdom));

		return sfp_ptr;
	}

	void output_linearization_domain(const positional_vdomain& dom,
			continuous::support_function_provider::const_ptr& X0,
			continuous::support_function_provider::const_ptr& X_next,
			math::lin_constraint_system<scalar_type>& domain_constraints) const {
		using namespace continuous;
		static int counter = 0;
		std::ofstream my_file;
		std::string filename_prefix = "/tmp/out_nonlin_" + to_string(counter);
		std::string filename_ini = filename_prefix + "_ini.txt";
		std::string filename_next = filename_prefix + "_next.txt";
		std::string filename_dom = filename_prefix + "_dom.txt";
		std::list<variable_id> v_list;
		v_list.push_back(dom.get_variable(0).get_id());
		v_list.push_back(dom.get_variable(1).get_id());
		my_file.open(filename_ini.c_str());
		{
			io::GEN_formatter formatter(my_file, 0);
			formatter.set_output_variables(v_list);
			formatter.output(*X0);
		}
		my_file.close();
		my_file.open(filename_next.c_str());
		{
			io::GEN_formatter formatter(my_file, 0);
			formatter.set_output_variables(v_list);
			formatter.output(*X_next);
		}
		my_file.close();
		my_file.open(filename_dom.c_str());
		{
			io::GEN_formatter formatter(my_file, 0);
			formatter.set_output_variables(v_list);
			constr_polyhedron<scalar_type> poly;
			poly.add_constraints(domain_constraints);
			formatter.output(poly);
		}
		my_file.close();
		IFLOGGER(DEBUG7) {
			int ret = system(
					std::string(
							"graph -TX -C -B -m2 -q-1 " + filename_dom
									+ " -s -m3 -q0.5 " + filename_ini
									+ " -s -m1 -q-1 " + filename_next).c_str());
		}
		++counter;
	}

	/** Obtains a polyhedron (in space-time) from a spacetime_flowpipe
	 *
	 * The quantification over the time variable and the last discrete
	 * jump are taken care of by mapping. */
	void get_poly_from_flowpipe(
			continuous::continuous_set_const_ptr& cset) const {
		using namespace continuous;
		using namespace support_function;

		// -------------------------------------------------------------
		// if cset is a flowpipe piece, compute its outer polytope and map it
		// -------------------------------------------------------------
		if (typename flowpipe::const_ptr cflow = boost::dynamic_pointer_cast<
				const flowpipe>(cset)) {
			// get a single outer polytope
			poly_ptr poly = cflow->compute_convex_outer_polyhedron(my_param.affine_error);
			// lock outer approximation so that containment checking can be performed
			flowpipe* nonconst_cflow = const_cast<flowpipe*>(cflow.get());
			nonconst_cflow->set_refinement_locked();

			// Now map the polyhedra (in space-time)
			if (!math::definitely(poly->is_empty())) {
//std::cout << std::endl << "Poly: " << **polys.begin() << std::endl;
				support_function_provider::const_ptr mapped_set;
				if (cflow->get_map()) {
//std::cout << "Map: " << *cflow->get_map() << std::endl;
					mapped_set = map_and_bloat(poly, *cflow->get_map(),
							cflow->get_bloating());
				} else {
					mapped_set = bloat<scalar_type>(poly,
							cflow->get_bloating());
				}
				IFLOGGER(DEBUG5) {
					hyperbox<scalar_type> bb =
							compute_bounding_box<scalar_type>(*mapped_set);
					LOGGER_OS(DEBUG5, __FUNCTION__) << "bounding box: " << bb
							<< std::endl;
				}
				cset = mapped_set;
			} else {
				cset = continuous::continuous_set_const_ptr();
			}
		}
	}

	/** Helper function to intersect a support function provider with a hyperbox */
	static continuous::support_function_provider::const_ptr intersect(
			const continuous::support_function_provider::const_ptr& U,
			const continuous::finite_hyperbox<scalar_type>& dbounds) {
		using namespace continuous;
		using namespace support_function;
		if (const continuous::finite_hyperbox<scalar_type>* h = dynamic_cast<const continuous::finite_hyperbox<scalar_type>*>(U.get())) {
			typename continuous::finite_hyperbox<scalar_type>::ptr b(
					new continuous::finite_hyperbox<scalar_type>());
			*b = compute_intersection(*h, dbounds);
			// @todo Check why this gives empty if polyhedron intersection is used
			return b;
		} else {
			throw std::runtime_error(std::string(__FUNCTION__)+": missing implementation");
		}
		return U;
	}

	/** Helper function to intersect a support function provider with a hyperbox */
	static continuous::support_function_provider::const_ptr intersect(
			const continuous::support_function_provider::const_ptr& U,
			const math::lin_constraint<scalar_type>& con) {
		using namespace continuous;
		using namespace support_function;
		if (const continuous::polyhedron<scalar_type>* h = dynamic_cast<const continuous::polyhedron<scalar_type>*>(U.get())) {
			typename continuous::polyhedron<scalar_type>::ptr h_clone(h->clone());
			h_clone->add_constraint(con);
			// @todo Check why this gives empty if polyhedron intersection is used
			return h_clone;
		} else {
			throw std::runtime_error(std::string(__FUNCTION__)+": missing implementation");
		}
		return U;
	}

	struct waiting_item {
		typedef continuous::support_function_provider::const_ptr set_const_ptr;
		set_const_ptr X;
		double time_offset;
		double previous_dwell;
	};

	struct waiting_comp {
	public:
		/** Returns true if lhs of less priority */
		bool operator() (const waiting_item& lhs, const waiting_item&rhs) const {
			return lhs.time_offset > rhs.time_offset;
		}
	};

	/** Compute time elapse from initial states X0 for nonlinear dynamics dp_nonlin
	 * under invariant constraints inv_constraints.
	 *
	 * Returns true if all successors have been exhausted
	 */
	bool time_elapse_nonlinear(
			continuous::continuous_set_collection& csets,
			continuous::support_function_provider::const_ptr X0,
			const typename continuous::typed_dynamics<scalar_type>& dyn,
			const math::lin_constraint_system<scalar_type>& inv_constraints, double time_offset = 0.0) const {
		using namespace continuous;
		using namespace support_function;

		bool all_exhausted = false;

		double time_horizon = get_time_horizon();
		double previous_dwell = time_horizon;

		typedef continuous::support_function_provider::const_ptr set_const_ptr;
		std::priority_queue<waiting_item,std::vector<waiting_item>,waiting_comp> waiting;
		waiting_item item;
		item.X = X0;
		item.time_offset = time_offset;
		item.previous_dwell = time_horizon;
		waiting.push(item);

		affine_map M;
		const positional_vdomain& dom = dyn.codom();
		// repeat until invariant violated or time horizon exhausted
		//support_function_provider::const_ptr X0 = support_function_provider::const_ptr(csup->clone());
		while (!all_exhausted) {
			bool exhausted = false;
			item = waiting.top();
			waiting.pop();
			X0 = item.X;
			time_offset = item.time_offset;
			previous_dwell = item.previous_dwell;

			LOGGER_OS(DEBUG5, __FUNCTION__) << "linearization starting from t="
					<< time_offset;

			continuous::support_function_provider::const_ptr U;

			// compute domain dom_inv, dynamics M, inputs U
			math::lin_constraint_system<scalar_type> domain_constraints;
			double nl_precision = my_param.linearization_error_factor * my_param.affine_error.abs();

			// error
			typename flowpipe::error_type fp_err = my_param.affine_error;

			affine_interpolation_result::type res;
			bool accepted;
			do {
				accepted = true;
				LOGGER_OS(DEBUG5, __FUNCTION__) << "interpolation with error " << nl_precision << " and affine flowpipe error " << fp_err;
				res = affine_interpolation<
					scalar_type>(dyn, *X0, M, U, domain_constraints,
					nl_precision);
				if (res != affine_interpolation_result::SUCCESS) {
					accepted = false;
				}

				// measure distance to the domain_constraints
				if (accepted) {
					scalar_type distance = distance_to_constraints(*X0,
							domain_constraints);
					scalar_type dstance_tolerance = my_param.domain_distance_factor
											* fp_err.abs();
					if (distance < dstance_tolerance) {
						LOGGER_OS(DEBUG6, __FUNCTION__) << "distance to domain "
								<< distance << " smaller than tolerance "
								<< dstance_tolerance << " for precision "
								<< nl_precision;
						accepted = false;
					}
				}
				if (!accepted) {
					if (my_param.split_linearization) {
						// get the center of X0
						finite_hyperbox<scalar_type> Xbounds =
								finite_bounding_box<scalar_type>(*X0);
						// compute the bounds on the error over the initial set
						finite_hyperbox<scalar_type> dbounds =
								derivative_err_bounds<scalar_type>(dyn, M, *X0);
						// split X0 along the largest error
//						const math::vdom_vector<scalar_type>& err_vec =
//								dbounds.get_g_dom();
						const math::vdom_vector<scalar_type>& err_vec =
								Xbounds.get_g_dom();
						unsigned int argmax_index;
						math::max(err_vec.get_vector(), argmax_index);
						// construct two constraints, one LE and one GE
						math::vdom_vector<scalar_type> pos_dir(dom);
						pos_dir[argmax_index] = scalar_type(1); // direction x_i
						scalar_type c_i = Xbounds.get_c_dom()[argmax_index]; // center in coordinate x_i -> c_i
						// constraint x_i >= c_i => -x_i + c_i <= 0
						math::lin_constraint<scalar_type> ge_con(-pos_dir, c_i,
								LE);
						// constraint x_i <= c_i =>  x_i - c_i <= 0
						math::lin_constraint<scalar_type> le_con(pos_dir, -c_i,
								LE);

						// intersect X0 with both constraints
						set_const_ptr X0_ge = intersect(X0, ge_con);
						set_const_ptr X0_le = intersect(X0, le_con);

						LOGGER_OS(DEBUG7, __FUNCTION__) << "splitting X with "
								<< ge_con << " and " << le_con;

						// continue with X0_ge, add X0_le to waiting list
						X0 = X0_ge;
						item.X = X0;
						waiting_item other = item;
						other.X = X0_le;
						waiting.push(other);
					} else {
						nl_precision = nl_precision
								* my_param.linearization_upscale_factor;
						fp_err = fp_err
								* sqrt(my_param.linearization_upscale_factor);
					}
				}
			} while (!accepted);

			if (res == affine_interpolation_result::SUCCESS) {

				/** Refine the bounds on U based on measurements (affine arithmetic) */
				if (my_param.refine_linearization_error) {
					// compute the bounds on the error over the whole poly
					constr_polyhedron<scalar_type> dom_poly;
					dom_poly.add_constraints(domain_constraints);

					finite_hyperbox<scalar_type> dbounds = derivative_err_bounds<
							scalar_type>(dyn,M, dom_poly);
					LOGGER_OS(DEBUG6,__FUNCTION__)
							<< "deriv err bounds in domain: " << dbounds << " versus " << nl_precision;
					U = intersect(U,dbounds);
					LOGGER_OS(DEBUG6,__FUNCTION__)
							<< "new err bounds in domain: " << U;
				}

				constr_polyhedron<scalar_type> dom_poly;
				dom_poly.add_constraints(domain_constraints);

				double estimated_dwell = estimate_max_dwell_time(ode_affine_dynamics<scalar_type>(M),U,*X0,dom_poly);
				if (estimated_dwell < 0.0)
					estimated_dwell = get_time_horizon()*my_param.default_dwell_time_factor;
				// estimated_dwell = std::min(estimated_dwell,get_time_horizon()/10);
				// this seems to be always worse: double nl_estimated_dwell = estimate_max_dwell_time(dyn,U,*X0,dom_poly);
				time_horizon = std::min(get_time_horizon(),time_offset+estimated_dwell);

				LOGGER_OS(DEBUG4,__FUNCTION__) << "time horizon " << time_horizon << ", max flow duration " << time_horizon-time_offset;

//									LOGGER_OS(DEBUG6,__FUNCTION__)
//											<< "dynamics " << *dp_nonlin
//											<< " for ini states " << *X0
//											<< " and precision " << nl_precision
//											<< std::endl;
//									LOGGER_OS(DEBUG6,__FUNCTION__)
//											<< " got affine approx dynamics "
//											<< M << " with input set " << U
//											<< " and domain of validity "
//											<< domain_constraints << std::endl;

				// instantiate a flowpipe object
				typename flowpipe::ptr sfp_ptr = create_flowpipe_object(X0,
						inv_constraints, M, U, time_offset,time_horizon);

				double t_last;
				if (!my_param.nonlinear_flowcut) {
					// restrict it to the domain
					restrict_to_sat_domain(sfp_ptr, domain_constraints,fp_err);
					t_last = sfp_ptr->get_time_domain().upper();

					// repeat with refined U
					if (my_param.refine_linearization_with_reach) {
						// for comparison
						finite_hyperbox<scalar_type> domain_box=finite_bounding_box<scalar_type>(dom_poly);
						// compute bounding box of reachable
						finite_hyperbox<scalar_type> reach_box = compute_finite_bounding_box<scalar_type>(*sfp_ptr);
						LOGGER_OS(DEBUG6,__FUNCTION__)
								<< "reach box: " << reach_box
								<< " versus " << domain_box;

						// compute error restricted to reachable box
						finite_hyperbox<scalar_type> dbounds =
								derivative_err_bounds<scalar_type>(dyn, M,
										reach_box);
						LOGGER_OS(DEBUG6,__FUNCTION__)
								<< "refined deriv err bounds in domain: " << dbounds
								<< " versus " << nl_precision;
						U = intersect(U, dbounds);
						LOGGER_OS(DEBUG6,__FUNCTION__)
								<< "refined new err bounds in domain: " << U;

						typename flowpipe::ptr sfp_ptr = create_flowpipe_object(X0,
								inv_constraints, M, U, time_offset,time_horizon);
						restrict_to_sat_domain(sfp_ptr, domain_constraints,fp_err);
					}
				} else {
					// restrict it to the domain
					restrict_to_invariant(sfp_ptr, domain_constraints,fp_err);
				}

				// Must compute t_last before intersecting with the invariant
				// so that we know when we've exhausted the time horizon even
				// though it cannot be exhausted inside invariant

				t_last = sfp_ptr->get_time_domain().upper();
				// If we've reached the time horizon, we're done
				if (math::maybe(
						math::numeric::is_GE(t_last, get_time_horizon()))) {
					exhausted = true;
				}
				// Check if progress has been made
				if (math::maybe(math::numeric::is_LE(t_last, time_offset))) {
					basic_warning(__FUNCTION__,
							"could not make time progress in time elapse of nonlinear dynamics",
							basic_warning::INCOMPLETE_OUTPUT);
					exhausted = true;
				}

				// restrict the flowpipe to the invariant
				const positional_vdomain& dom = M.codomain();
				math::lin_constraint_system<scalar_type> state_invariant_cons =
						compute_state_variable_constraints<scalar_type>(dom,
								inv_constraints);
				restrict_to_invariant(sfp_ptr, state_invariant_cons,fp_err);

				// if invariant reduce time domain, we're done
				// @todo: stop when invariant is violated by all states
//				if (math::definitely(
//						math::numeric::is_LT(sfp_ptr->get_time_domain().upper(),
//								t_last))) {
//					exhausted = true;
//				}

				add_user_template_directions(sfp_ptr);

				// add flowpipe piece to result
				csets.insert(sfp_ptr, false); // false <- no redundancy check

				if (!my_param.nonlinear_flowcut) {
					// take the last states inside the domain
					support_function_provider::const_ptr X_next =
							sfp_ptr->get_last_states();
					t_last = sfp_ptr->get_time_domain().upper();

					IFLOGGER(DEBUG6) {
						output_linearization_domain(dom, X0, X_next,
								domain_constraints);
					}

					// prepare next iteration
					if (!X_next->is_empty() && !exhausted) {
						item.X = X_next;
						item.time_offset = t_last;
						item.previous_dwell = t_last - time_offset;
						waiting.push(item);
					}
				} else {
					// intersect at the boundary of the linearization domain
					previous_dwell = -1; // find the min dwell

					// for each constraint of the boundary
					for (typename math::lin_constraint_system<scalar_type>::const_iterator it =
							domain_constraints.begin();
							it != domain_constraints.end(); ++it) {
						// flip the constraint
						math::lin_constraint<scalar_type> con(
								-it->get_canonic_l(), it->get_canonic_sign());
						LOGGER_OS(DEBUG6, __FUNCTION__)
								<< "intersecting with domain boundary " << con;
						// clone the flowpipe
						typename flowpipe::ptr sfp_copy(sfp_ptr->clone());

						// intersect the constraint
						typedef typename flowpipe::time_interval time_interval;
						typedef std::vector<time_interval> intv_list_type;
						intv_list_type intv_list =
								sfp_copy->get_maybe_satisfying_subdomains(con,
										my_param.affine_error,false);
						if (!intv_list.empty()) {
							time_interval intersect_intv = *intv_list.begin();
							LOGGER_OS(DEBUG6, __FUNCTION__)
									<< "found time interval "
									<< intersect_intv.lower() << ","
									<< intersect_intv.upper();

							sfp_copy->restrict_to_subdomain(intersect_intv);
							t_last = intersect_intv.lower();

							//pick the smallest delta as previous_dwell
							double delta = t_last - time_offset;
							if (previous_dwell < 0.0 || delta < previous_dwell) {
								previous_dwell = delta;
							}

							// restrict_evolutions = true
							sfp_copy->add_outer_constraint(con, true);
							// get a single outer polytope
							typename flowpipe::cut_point_method m;
							m.type =
									flowpipe::cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES;
							m.piece_count = 1;
							m.approx_error = my_param.affine_error;
							polyhedron_collection<scalar_type> polys =
									sfp_copy->compute_outer_polyhedra(m);

							// Now map the polyhedra (in space-time)
							typename polyhedron_collection<scalar_type>::element_type poly =
									*polys.begin();
							if (!math::definitely(poly->is_empty())) {
								size_t orig_nb =
										poly->get_constraints()->size();
								poly->remove_redundant_constraints();
								size_t before_proj_nb =
										poly->get_constraints()->size();
								variable_id_set vis;
								vis.insert(
										sfp_copy->get_time_variable().get_id());
								poly->existentially_quantify_variables(vis);
								size_t after_proj_nb =
										poly->get_constraints()->size();
								poly->remove_redundant_constraints();
								size_t final_nb =
										poly->get_constraints()->size();
								LOGGER_OS(DEBUG6, __FUNCTION__) << "poly orig "
										<< orig_nb << " nred " << before_proj_nb
										<< " proj " << after_proj_nb
										<< " final " << final_nb << ": "
										<< poly;

								support_function_provider::const_ptr X_next =
										poly;

								// csets->insert(poly, false);
								IFLOGGER(DEBUG6) {
									output_linearization_domain(dom, X0, X_next,
											domain_constraints);
								}

								// prepare next iteration
								if (math::definitely(
										math::numeric::is_LT(t_last, get_time_horizon()))) {
									item.X = X_next;
									item.time_offset = t_last;
									item.previous_dwell = t_last - time_offset;
									waiting.push(item);
								}
							}
						}
					}
				}
				LOGGER_OS(DEBUG4, __FUNCTION__) << "covered " << previous_dwell << " time units, up to time "
						<< waiting.top().time_offset << ", " << waiting.size() << " waiting";
				all_exhausted = waiting.empty();
			} else {
				// requires splitting
				basic_warning(__FUNCTION__,
						"could not make time progress in time elapse of nonlinear dynamics",
						basic_warning::INCOMPLETE_OUTPUT);
				all_exhausted = true;
			}
		}

		return all_exhausted;
	}


	/** Compute time elapse from initial states X0 for nonlinear dynamics dp_nonlin
	 * under invariant constraints inv_constraints.
	 */
	continuous::continuous_set_ptr time_elapse_nonlinear(
			continuous::support_function_provider::const_ptr X0,
			const typename continuous::typed_dynamics<scalar_type>& dyn,
			const math::lin_constraint_system<scalar_type>& inv_constraints) const {
		using namespace continuous;
		using namespace support_function;

		LOGGER(DEBUG4, "continuous_post_stc::post",
				"time elapse for nonlinear dynamics");

		double time_offset = 0.0;
		continuous_set_collection::ptr csets(new continuous_set_collection());
		bool exhausted = time_elapse_nonlinear(*csets,X0,dyn,inv_constraints,time_offset);

		return csets;
	}

	/** Compute time elapse for affine dynamics */
	continuous::continuous_set_ptr time_elapse_affine(
			continuous::support_function_provider::const_ptr& csup,
			const continuous::ode_affine_dynamics<scalar_type>& M,
			typename continuous::polyhedron<scalar_type>::const_ptr new_inv, double time_horizon) const {
		using namespace continuous;
		using namespace support_function;

		LOGGER(DEBUG4, "continuous_post_stc::post",
				"time elapse for affine dynamics");

		// divide the dynamics into state equation and additive input disturbance,
		// given the invariant (over both state and input sets)
		typename polyhedron<scalar_type>::const_ptr red_inv = new_inv;
//		ode_affine_dynamics<scalar_type> M = convert_to_ODE(dyn,new_inv,red_inv);

		// affine dynamics, produce a single set
		support_function_provider::const_ptr U = M.get_U();
		math::lin_constraint_system<scalar_type> red_inv_constraints;
		if (red_inv) {
			red_inv_constraints = *red_inv->get_constraints();
		}

		IFLOGGER(DEBUG6) {
			LOGGER_OS(DEBUG6,__FUNCTION__) << "dynamics: " << M << ", invariant: " << red_inv_constraints;
			// report eigenvalues
			IFLOGGER(DEBUG6)
			{
				using namespace math;
				vector<scalar_type> v_real_d, v_imag_d;
				matrix<scalar_type> V_d;
				compute_eigenvalues(M.get_A().get_matrix(), v_real_d, v_imag_d, V_d);
				LOGGER_OS(DEBUG6, __FUNCTION__) << "real parts of Eigenvalues:"
						<< v_real_d << ", imag parts of Eigenvalues:" << v_imag_d;
				//LOGGER_OS(DEBUG4, __FUNCTION__) << " Eigenvectors:" << V_d;
			}
		}

		// create the flowpipe object
		double time_offset = 0.0;
		typename flowpipe::ptr sfp_ptr = create_flowpipe_object(csup,
				red_inv_constraints, M, U, time_offset, time_horizon);

		// restrict the flowpipe to the invariant
		restrict_to_invariant(sfp_ptr, red_inv_constraints,my_param.affine_error);
		add_user_template_directions(sfp_ptr);
		return sfp_ptr;
	}

	/** Return the continuous_set that results from applying time elapse
	 * with the time constraints tcons to the continuous set *cset. */
	virtual continuous::continuous_set_ptr post(const time_constraints& tcons,
			const continuous::continuous_set_const_ptr& cset_orig) const {
		LOGGERSW(DEBUG1, "continuous_post_stc::post",
				"Continuous post with continuous_post_stc");
		using namespace continuous;
		using namespace support_function;
		using namespace math;
		using namespace math::numeric;

		continuous::continuous_set_const_ptr cset(cset_orig);
		typedef ode_affine_dynamics<scalar_type> ode_aff_dyn;
		continuous_set_ptr ret_set = continuous_set_ptr();
		// -------------------------------------------------------------
		// if cset is a flowpipe piece, compute its outer polytope and map it
		// -------------------------------------------------------------
		get_poly_from_flowpipe(cset);
		if (cset)
			if (support_function_provider::const_ptr csup =
					boost::dynamic_pointer_cast<const support_function_provider>(
							cset)) {

				typename polyhedron<scalar_type>::const_ptr inv_poly_ptr;
				if (tcons.get_invariant()) {
					inv_poly_ptr = boost::dynamic_pointer_cast<
							const polyhedron<scalar_type> >(
							tcons.get_invariant());
					if (!inv_poly_ptr) {
						std::string tname =
								typeid(*tcons.get_invariant().get()).name();
						throw basic_exception(
								"continuous_post_stc::post cannot handle invariant of type "
										+ tname);
					}
				}

				typename ode_aff_dyn::const_ptr dp =
						boost::dynamic_pointer_cast<const ode_aff_dyn>(
								tcons.get_dynamics());
				typename typed_dynamics<scalar_type>::const_ptr dp_nonlin =
						boost::dynamic_pointer_cast<
								const typed_dynamics<scalar_type> >(
								tcons.get_dynamics());
				if (dp || dp_nonlin) {
					LOGGER(DEBUG3, "continuous_post_stc::post",
							"computing support_function_time_elapse");

					// first, we treat degenerated cases
					if ((dp && dp->is_universe())) {
						// there is no restriction on the dynamics,
						// so the entire invariant is reachable.
						// simply return a copy of the invariant

						LOGGER(DEBUG3, "continuous_post_stc::post",
								"no flow restriction (flow predicate is always true)");
						ret_set = continuous_set_ptr(inv_poly_ptr->clone());
					} else if ((dp && dp->is_zero())) {
						// there is no change in the variables,
						// since the derivatives are identical zero.
						// Simply return a copy of the initial states

						LOGGER(DEBUG3, "continuous_post_stc::post",
								"no time elapse in location (unsat flow predicate)");
						ret_set = continuous_set_ptr(cset_orig->clone());
					} else {
						// call time elapse calculation
						// restrict invariant with const variables
						typename polyhedron<scalar_type>::const_ptr new_inv =
								inv_poly_ptr;
						if (false && new_inv && dp) {
							const affine_map& M_orig(*dp);
							// First, construct the bounding box of the derivatives in this invariant
							support_function::sf_unary<scalar_type> deriv_set(
									new_inv, M_orig);
							// Go through the variables in M_orig and find the ones with derivative zero
							variable_id_set constvars = get_vars_bound_to_zero<
									scalar_type>(deriv_set);

							if (!constvars.empty()) {
								std::stringstream ss;
								logger::copyfmt_to(ss);
								print_variable_id_set(ss, constvars);
								LOGGER(DEBUG4, "continuous_post_stc",
										"adding const constraints to variables "+ss.str());
								//std::cout << "constant state vars:"; print_variable_id_set(std::cout,constvars); std::cout << std::endl;
								// Get bounds on the initial values of constvars
								positional_vdomain const_dom(constvars);
								typedef vector<scalar_type> vector_type;
								std::list<vector_type> constdirections;
								support_function::choose_directions(const_dom,
										constdirections); // R is the number of directions.
								typedef std::set<vdom_vector<scalar_type>,
										numeric::lex_comp_less<
												scalar_type, vdom_vector> > dir_set_type;

								dir_set_type dir_set;

								for (typename std::list<
										vector<scalar_type> >::const_iterator it =
										constdirections.begin();
										it != constdirections.end(); ++it) {
									vdom_vector<scalar_type> d(const_dom,
											*it);
									dir_set.insert(d);
								}

								typename constr_polyhedron<scalar_type>::ptr poly(
										new constr_polyhedron<scalar_type>(
												compute_outer_poly<scalar_type>(
														*csup, dir_set)));
								IFLOGGER(DEBUG5) {
									ss.str("");
									ss << poly;
									LOGGER(DEBUG5, "continuous_post_stc",
											"const constraints are: "+ss.str());
								}
								// add poly to new_inv
								poly->add_constraints(
										new_inv->get_constraints());
								new_inv = poly;
								//std::cout << "new invariant:" << new_inv;
							}
						}

						if (dp) {
							// if there are extra variables csup, project
							variable_id_set vars = csup->get_variable_ids();
							//variable_id_set target_nonstate_variables = inst.tinv_ptr->get_variable_ids();
							variable_id_set state_variables = dp->domain().get_variable_ids();
							set_difference_assign(vars, state_variables);

							if (!vars.empty()) {
								LOGGER(DEBUG7, __FUNCTION__,
										"removing nonstate "+to_string(vars.size())+" variables for time elapse ");
								positional_vdomain cdom(csup->get_variable_ids());
								affine_map M_proj = affine_map::projection_map(dp->domain(),cdom);
								csup = support_function_provider::const_ptr(
										new sf_unary<scalar_type>(csup, M_proj));
							}

							// check if nonzero time elapse is possible
							bool use_time_horizon_estimation = true;
							bool skip_time_elapse_if_horizon_zero = true;

							double time_horizon = get_time_horizon();
							if (use_time_horizon_estimation && new_inv) {
								double estimated_dwell = -1.0;
								const support_function_provider& test = *csup;
								estimated_dwell = estimate_max_dwell_time(*dp, dp->get_U(), *csup,
										*new_inv);
								if (!definitely(is_LT(estimated_dwell,0.0)) && definitely(is_LT(estimated_dwell,time_horizon))) {
									time_horizon = estimated_dwell;
									// don't go below zero
									time_horizon = std::max(time_horizon,0.0);
									LOGGER(DEBUG, __FUNCTION__,
											"reduced time horizon to " + to_string(time_horizon));
								}
							}

							bool compute_time_elapse = true;
							if (skip_time_elapse_if_horizon_zero && is_MEQ(time_horizon, 0.0)) {
								compute_time_elapse = false;

								// if it's a flowpipe check first whether it's to be
								// mapped and bloated
								if (typename flowpipe::const_ptr p = boost::dynamic_pointer_cast<const flowpipe>(csup)) {
									// If to be bloated, we need to use the outer poly
									if (p->get_bloating()) {
										compute_time_elapse = true;
									}
								}
							}
							if (compute_time_elapse){
								ret_set = time_elapse_affine(csup, *dp, new_inv,
										time_horizon);
							}
							if (!compute_time_elapse) {
								ret_set = continuous_set_ptr(cset_orig->clone());
							}
						} else if (dp_nonlin) {
							lin_constraint_system<scalar_type> inv_constraints;
							if (new_inv)
								inv_constraints.push_back(
										*new_inv->get_constraints());
							return time_elapse_nonlinear(csup, *dp_nonlin,
									inv_constraints);
						}
					}
				} else {
					//				std::cerr << "offended by dynamics:"
					//						<< tcons.get_dynamics()->get_predicate() << std::endl;
					std::string tname =
							typeid(*tcons.get_dynamics().get()).name();
					throw std::runtime_error(
							"continuous_post_stc::post cannot handle dynamics of type "
									+ tname);
				}
			} else {
				std::string tname = typeid(*cset.get()).name();
				throw std::runtime_error(
						"continuous_post_stc::post cannot handle continuous_set of type"
								+ tname);
			}

		return ret_set;
	}
	;

private:
	time_elapse_parameters my_param;

};

}

#endif /* CONTINUOUS_POST_STC_H_ */
