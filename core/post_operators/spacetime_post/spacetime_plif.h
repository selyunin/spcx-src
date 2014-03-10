/*
 * spacetime_plif.h
 *
 *  Created on: Oct 18, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_PLIF_H_
#define SPACETIME_PLIF_H_

#include "math/numeric/rel_abs_error.h"
#include "math/numeric/interval.h"
#include "math/numeric/basic_geometry.h"
#include "math/vdom/lin_constraint_system.h"
#include "math/unique_vector_to_value_store.h"
#include "math/unique_scalar_to_value_store.h"

#include "core/continuous/continuous_dynamics/ode_affine_dynamics.h"
#include "core/continuous/support_function_provider.h"
#include "core/continuous/polyhedra/hyperbox/finite_hyperbox_utility.h"
#include "core/continuous/support_function/sf_base/sf_unary_ref.h"

#include "core/continuous/polyhedra/polyhedron_collection.h"
#include "math/vdom/positional_vdomain.h"

#include "extern/plif/plif.h"

namespace spacetime {

using namespace continuous;
using namespace math;
using namespace math::numeric;

/** A structure for representing affine dynamics with initial and boundary conditions
 *
 * The dynamics are
 *  \f$ \dot x = Ax + b + u\f$, where \f$ u \in U \subseteq R^n\f$ is a set of nondeterministic
 *  inputs. The initial conditions are \f$ x \in X_0\f$ and the boundary conditions are
 *  \f$ x(t) \in V \f$.
 *
 *  Null pointers for dyn.get_U() and X0 are interpreted as the set consisting of a zero vector.
 * */
template<typename scalar_type>
struct affine_support_flowpipe_description {
	ode_affine_dynamics<scalar_type> dyn;
	typename lin_constraint_system<scalar_type>::const_ptr V;
	support_function_provider::const_ptr X0;
};

/** A class for approximating flowpipes with piecewise linear support functions over time
 *
 *
 * The class is capable of computing support evolutions on demand up to a certain error bound,
 * retains already computed evolutions, and can refine computed evolutions up to a given error bound.
 */
class spacetime_plif {
public:
	typedef double scalar_type;
	typedef math::numeric::rel_abs_error<scalar_type> error_type;
	typedef scalar_type duration;
	//	typedef interval<duration> time_interval;
	//	typedef interval<scalar_type> scalar_interval;
	typedef plif::interval time_interval;
	typedef plif::interval scalar_interval;
	typedef plif::piecewise_linear_function plf_type;
	typedef plif::piecewise_linear_interval_function plif_type;
	typedef math::vdom_matrix<scalar_type> matrix_type;
	typedef math::vdom_vector<scalar_type> vector_type;
	typedef vector_type direction;
	typedef std::pair<plif_type, error_type> annotated_plif;
	typedef unique_vector_to_value_store<scalar_type, math::vdom_vector,
			annotated_plif> evolution_cache; // stores the evolutions for each direction
	typedef std::set<duration,math::numeric::comp_less<scalar_type> > time_point_set;

	/** Construct a flowpipe approximation given an affine flowpipe */
	spacetime_plif(affine_support_flowpipe_description<scalar_type> prob,
			time_interval intv);

	/** Returns the domain of the state variables of the flowpipe */
	const positional_vdomain& domain() const;

	/** Get the time domain over which all evolutions are defined */
	const time_interval& get_time_domain() const;

	/** Get or compute support evolution for given direction and error bound
	 *
	 * The evolution is cached for later retrieval. The direction is normed according to the norm defined by
	 * get_norm(). The evolution for an arbitrary (nonzero) direction d is therefore obtained with
	 * evo=get_or_compute_evolution(d,err)*get_norm(d);
	 * */
	const annotated_plif& get_or_compute_evolution(direction d,
			const error_type& err);

	/** Get a support evolution for given direction if already in cache, and null otherwise
	 *
	 * */
	const annotated_plif* get_evolution(direction d);

	/** Get the subdomains in which the constraint is maybe satisfied for all t
	 *
	 * A constraint is "maybe satisfied" by an interval if there is at least one point in the
	 * interval that satisfies the constraint.
	 * */
	std::vector<time_interval> get_maybe_satisfying_subdomains(const lin_constraint<scalar_type>& con, const error_type& err, bool widen_by_error);

	/** Get the subdomains in which the constraint is definitely satisfied for all t
	 *
	 * A constraint is "definitely satisfied" by an interval if all points in the
	 * interval satisfy the constraint.
	 * */
	std::vector<time_interval> get_definitely_satisfying_subdomains(const lin_constraint<scalar_type>& con, const error_type& err);

	/** Restrict the domain of *this (including cached evolutions) to a subdomain
	 * */
	void restrict_to_subdomain(const time_interval& intv);

	/** Restrict the evolutions to the support of the constraint system */
	void restrict_to_constraints(const lin_constraint_system<scalar_type>& con);

	/** Returns whether the outer approximation contains the flowpipe f
	 *
	 * If refine is true, then f is refined until the property can be decided.
	 */
	math::tribool decide_outer_contains(spacetime_plif& f, bool refine = false) const;

	/** The different methods for obtaining cut points
	 *
	 * If type=FIXED_COUNT_SAME_SIZE_PIECES, the number of pieces is specified
	 * by piece_count and all pieces are chosen as the same size.
	 * If type=TIME_STEP, the time step is specified by time_step.
	 * If type=MIN_CONCAVE_PIECES, the smallest number of pieces for
	 * given directional approximation error approx_error is computed.
	 * If type=ALL_PIECES, all concave pieces are returned. This is the most
	 * exact solution possible. */
	struct cut_point_method {
		typedef enum {
			FIXED_COUNT_SAME_SIZE_PIECES,
			TIME_STEP,
			MIN_CONCAVE_PIECES,
			ALL_PIECES
		} method_type;

		cut_point_method(method_type t = ALL_PIECES, size_t pc = 0,
				duration ts = (0.0), error_type aprerr = error_type(0.0, 0.0),
				bool lowisref = true) :
				type(t), piece_count(pc), time_step(ts), approx_error(aprerr), lower_is_error_reference(
						lowisref) {
		}
		;

		method_type type;
		size_t piece_count;
		duration time_step;
		error_type approx_error;
		bool lower_is_error_reference;
	};

	/** Replace the upper bounds in the evolution cache by their convex hull */
	void assign_convex_hull();

	/** Obtain a set of convex flowpipe that cover *this */
	std::vector<spacetime_plif> convexify(const cut_point_method& m) const;

	/** Get outer polyhedral approximation using a given method for dividing into pieces
	 *
	 * If no evolutions have yet been computed, a universe polyhedron is returned. */
	polyhedron_collection<scalar_type> compute_outer_polyhedra(
			const cut_point_method& m) const;

	/** Return the number of evolutions stored in the cache */
	size_t size() const;

	/**	Set cut point computation method */
	void set_cut_point_method(cut_point_method m);

	/**	Get current cut point computation method */
	cut_point_method get_cut_point_method() const;

	/** Get the time variable of the set */
	const variable& get_time_variable() const;

	/** Obtain the norm of a direction */
	static scalar_type get_norm(const direction& d);

	/** Obtain the normed direction */
	static direction get_normed_direction(const direction& d);

	/** Returns the system description */
	const affine_support_flowpipe_description<scalar_type>& get_affine_support_flowpipe_description() const;

	/** Activates simplification of concave pieces in evolution computations */
	void set_simplify_concave(bool b);

	/** Activates simplification of convex pieces in evolution computations */
	void set_simplify_convex(bool b);

private:
	typedef finite_hyperbox<scalar_type> box_type; // for auxiliary sets

	/** A structure for all auxiliary variables depending on delta */
	struct delta_parameters {
		matrix_type Phi; //=e^(A*delta)
		matrix_type Phi_T; //=Phi^T
		matrix_type Phi1; //=A^{-1}(e^(A*delta)-I)
		matrix_type Phi2; //=A^{-2}(e^(A*delta)-I-delta*A)

		box_type E_X0f;
		box_type E_X0b;
		box_type E_U; // error for U without b
		support_function_provider::ptr negAPhi2U; // error for U without b
		box_type E_Ub; // error for U with b
		support_function_provider::ptr negAPhi2Ub; // error for U with b
	};
	typedef unique_scalar_to_value_store<scalar_type, delta_parameters>
			delta_parameters_cache_type;

	/** Compute the support function of S in direction l. */
	static scalar_type
	rho(const vector_type& l, const support_function_provider& S);

	/** Compute the support vector of S in direction l.
	 *
	 * Throws if the support vector cannot be computed by S. */
	vector_type
	rho_vec(const vector_type& l, const support_function_provider& S);

	/** Compute support evolution for normed direction, error bound, and initial time step
	 * */
	annotated_plif compute_evolution(const direction& d, const error_type& err);

	/** A structure for capturing an interval associated with a direction over a time interval
	 */
	struct timed_directional_interval {
		duration t; // time point in absolute time
		duration delta; // time duration
		direction d; // direction at time t
		scalar_interval itv; // interval
	};

	/** Comparing timed_directional_interval on time */
	struct compare_timed_directional_interval_times {
		bool operator()(const timed_directional_interval& a,
				const timed_directional_interval& b) {
			return definitely(is_LT(a.t, b.t));
		}
		;
	};

	/** A list of timed directional intervals
	 *
	 * This is using a simple list to keep overhead low when merging. */
	class timed_directional_interval_sequence {
	public:
		typedef std::list<timed_directional_interval>::const_iterator
				const_iterator;
		typedef std::list<timed_directional_interval>::iterator iterator;
		typedef std::list<timed_directional_interval>::reverse_iterator reverse_iterator;


		/** Construct an empty list */
		timed_directional_interval_sequence() {
		}
		;

		/** Construct from a single interval */
		explicit timed_directional_interval_sequence(
				const timed_directional_interval& p) {
			my_list.push_back(p);
		}

		/** Return the size of the sequence */
		size_t size() const {
			return my_list.size();
		};

		/** Add an interval at the end of the list
		 *
		 * The interval must no come before the time at the end of *this. */
		void push_back(timed_directional_interval& p) {
			//assert(!compare_timed_directional_interval_times()(p,*my_list.rbegin()));

			my_list.push_back(p);
		}

		/** Add an interval at the beginning of the list
		 *
		 * The interval must no come before the time at the end of *this. */
		void push_front(timed_directional_interval& p) {
			//assert(!compare_timed_directional_interval_times()(p,*my_list.rbegin()));

			my_list.push_front(p);
		}

		/** Concatenate with a temporary sequence
		 *
		 * The temp sequence is emptied by the operation.
		 * The temp sequence must start at a time point not before *this. */
		void concatenate_with_temp(timed_directional_interval_sequence& temp) {
			//assert(!compare_timed_directional_interval_times()(*temp.begin(),*my_list.rbegin()));
			//@todo if the sequences overlap, reconcile the overlapping breakpoints
			my_list.splice(my_list.end(), temp.my_list);
		}

		/** Concatenate with a temporary sequence and simplify
		 *
		 * The temp sequence is emptied by the operation.
		 * The temp sequence must start at a time point not before *this. */
		void concatenate_with_temp_and_simplify(timed_directional_interval_sequence& temp, const error_type& err) {
			//assert(!compare_timed_directional_interval_times()(*temp.begin(),*my_list.rbegin()));
			//@todo if the sequences overlap, reconcile the overlapping breakpoints
			iterator last_old_it=my_list.end();
			if (my_list.size()>0)
				--last_old_it;
			my_list.splice(my_list.end(), temp.my_list);
			if (last_old_it!=my_list.end()) {
				size_t old_size = my_list.size();
				simplify_concave(last_old_it,my_list.end());
				LOGGER(DEBUG7, __FUNCTION__, "simplified sequence of "+to_string(old_size)+ " to "+to_string(my_list.size()));
			}
		}

		/** Simplify by removing concave points that lie within error bound err
		 *
		 * This assumes that all points already satisfy the given error bound.
		 * Points will be eliminated such that the bound isn't increased.
		 * */
		void simplify_concave(iterator it_a,iterator it_end) {
			if (it_a != it_end) {
				iterator it_b = it_a;
				++it_b;
				if (it_b != it_end) {
					iterator it_c = it_b;
					++it_c;
					while (it_c != it_end) {
						// if a+,b+,c+ is convex and a-,b-,c- is concave, remove b
						if (is_concave(*it_a,*it_b,*it_c)) {
							// add b's delta to a
							it_a->delta += it_b->delta;
							// remove b
							it_b = my_list.erase(it_b);
							it_c = it_b;
							++it_c;
						} else {
							// continue to next point
							++it_a;
							++it_b;
							++it_c;
						}
					}
				}
			}
		};

		/** Const iterator to the beginning of the list */
		const_iterator begin() const {
			return my_list.begin();
		}

		/** Const iterator past the end of the list */
		const_iterator end() const {
			return my_list.end();
		}

		/** Iterator to the beginning of the list */
		iterator begin() {
			return my_list.begin();
		}

		/** Iterator past the end of the list */
		iterator end() {
			return my_list.end();
		}

		/** Iterator to the beginning of the reversed list */
		reverse_iterator rbegin() {
			return my_list.rbegin();
		}

		/** Iterator past the end of the reversed list */
		reverse_iterator rend() {
			return my_list.rend();
		}

		/** Sort the list */
		void sort() {
			my_list.sort(compare_timed_directional_interval_times());
		}

		/** Output to stream
		 *
		 * The points are separated by std::endl, the values are separated by spaces:
		 * time delta d itv.lower() itv.upper()
		 * The format of the direction and the interval depends on the corresponding
		 * operator<<.
		 */
		void print(std::ostream& os = std::cout) {
			for (const_iterator it = begin(); it != end(); ++it) {
				os << it->t << " " << it->delta << " " << it->d << " "
						<< it->itv.lower() << " " << it->itv.upper()
						<< std::endl;
			}
		}

		/** Returns true if a,b,c is definitely concave
		 *
		 * The graph of the points is concave if the upper bounds are convex and the lower bounds are concave */
		static bool is_concave(const timed_directional_interval& a,
				const timed_directional_interval& b,
				const timed_directional_interval& c) {
			using namespace math::numeric;
			return maybe(is_ccw(a.t,a.itv.upper(), b.t, b.itv.upper(), c.t, c.itv.upper())
					&& is_ccw(c.t, c.itv.lower(), b.t, b.itv.lower(), a.t, a.itv.lower()));
		}

	private:
		//std::set<timed_directional_interval,compare_timed_directional_interval_times> my_list;
		std::list<timed_directional_interval> my_list;
	};

	/** Return the linear interpolation of the interval sequence */
	static plif_type linear_interpolation(
			const timed_directional_interval_sequence& seq);

	/** Return the stepwise interpolation of the interval sequence
	 *
	 * This is equivalent to a zero order hold.*/
	static plif_type left_step_interpolation(
			const timed_directional_interval_sequence& seq);

	/** Create a sequence of support values on Omega_X0 with given error
	 *
	 * The error at time t is below err_begin, at t+delta below err_end.
	 * The direction vector is only computed at the beginning of each Omega.
	 * */
	timed_directional_interval_sequence omega_X0_err(const direction& d,
			const duration& t, const duration& delta,
			const error_type& err_begin, const error_type& err_end);

	/** Create a sequence of support values on Omega_X0 with given error
	 *
	 * The error at time t is below err_begin, at t+delta below err_end.
	 * The direction vector is only computed at the beginning of each Omega.
	 * */
	struct omega_X0_err_state {
		vector_type d_end;
		vector_type vec_end;
		scalar_type rho_end;
		scalar_type t;
		scalar_type delta;
		error_type  err_end;
	};

	timed_directional_interval_sequence omega_X0_err_simplified(const direction& d,
			const duration& t, const duration& delta,
			const error_type& err_begin, const error_type& err_end);

	/** Create a sequence of time intervals with given error on Omega
	 *
	 * The error at time t is below err_begin, at t+delta below err_end.
	 * The direction vector is only computed at the beginning of each Omega.
	 * */
	timed_directional_interval_sequence omega_X0_delta(const direction& d,
			const duration& t, const duration& delta,
			const error_type& err_begin, const error_type& err_end, const vector_type& vec_Omega_beg,
			vector_type& vec_Omega_end);

	/** Create a sequence of support values on Omega_Ub with given error
	 *
	 * The error at time t is below err_begin, at t+delta below err_end.
	 * The direction vector is only computed at the beginning of each Omega.
	 *
	 * The resulting sequence refers to the integral over U at the *end* of the time interval
	 * -- the value at the beginning being zero.
	 * */
	timed_directional_interval_sequence omega_Ub_err(const direction& d,
			const duration& t, const duration& delta,
			const error_type& err_begin, const error_type& err_end);

	/** Compute error bounds for rho_Omega(delta) */
	static void compute_rho_Omega_err(scalar_type& err_low, scalar_type& err_upp,
			const vector_type& ef, const vector_type& eb, const vector_type& d,
			const duration& lambda, const scalar_type& err_quad_psi_upp,
			const scalar_type& err_lin_psi_low);

	/** Compute rho_{X0}(d)
	 */
	scalar_type rho_X0(const vector_type& d);

	/** Compute the support vector for rho_{X0}(d) */
	vector_type rho_vec_X0(const vector_type& d);

	/** Get the variable domain */
	const positional_vdomain& domain();

	/** Compute rho_{U+b}(d)
	 */
	scalar_type rho_Ub(const vector_type& d);

	/** Compute rho_{E_Omega(delta,t)}(d)
	 *
	 * lambda = t/delta
	 */
	static scalar_type
	rho_E_Omega(const vector_type& ef, const vector_type& eb,
			const vector_type& d, const duration& lambda);

	/** Return the value of psi at time t+delta up to a given error
	 *
	 *	The direction d is the direction at the beginning of the time interval I.
	 *  */
	scalar_interval psi_err_value(const direction& d, const duration& t,
			const duration& delta, const error_type& err_begin,
			const error_type& err_end);

	/** Create a sequence of interval samples of Psi for a given direction, in a given time interval up to a given error
	 *
	 *	The time interval is [t,t+delta].
	 *	The direction d is the direction at time t.
	 *  */
	timed_directional_interval_sequence psi_err(const direction& d,
			const duration& t, const duration& delta,
			const error_type& err_begin, const error_type& err_end);

	/** Create a sequence of time intervals with given error on Psi
	 *
	 * If use_rate is not specified or true, then the resulting error
	 * is below err_end-err_begin (To avoid infinite recursion, ensure err_end > err_begin.)
	 * Otherwise, the resulting error is below min(err_begin,err_end).
	 * (To avoid infinite recursion, ensure err_begin and err_end are strictly greater than
	 * zero.)
	 * */
	timed_directional_interval_sequence
	psi_delta(const direction& d, const duration& t, const duration& delta,
			const error_type& err_begin, const error_type& err_end,
			bool use_Ub, bool use_rate = true);

	/** Returns whether b or U are nonzero */
	bool has_Ub() const;

	/** Obtain E_X0f = E_Omega⁺(delta) */
	const box_type& E_X0f(const duration& delta) const;

	/** Obtain E_X0b = E_Omega⁻(delta) */
	const box_type& E_X0b(const duration& delta) const;

	/** Obtain E_U = E_Psi(U,delta) */
	const box_type& E_U(const duration& delta) const;

	/** Obtain negAPhi2U = -A*Phi2(A,delta)*U */
	const support_function_provider::ptr& negAPhi2U(const duration& delta) const;

	/** Obtain E_Ub = E_Psi(U+b,delta) */
	const box_type& E_Ub(const duration& delta) const;

	/** Obtain negAPhi2Ub = -A*Phi2(A,delta)*(U+b) */
	const support_function_provider::ptr
	& negAPhi2Ub(const duration& delta) const;

	/** Obtain e^{A\delta}^T */
	const matrix_type& exp_AdeltaT(const duration& delta) const;

	/** Obtain phi1(delta) */
	const matrix_type& phi1(const duration& delta) const;

	/** Returns -1 or the smallest time step that caused numerical problems */
	/**
	 * Find parameters for given delta
	 *
	 * To avoid looking up in the cache every time, check if the current pointer
	 * fits.
	 */
	const delta_parameters& get_delta_params(const duration& delta) const;

	/**
	 * Find parameters for given delta
	 *
	 * To avoid looking up in the cache every time, check if the current pointer
	 * fits.
	 */
	const delta_parameters& get_delta_params(const duration& delta);

	/** Compute the auxiliary sets for given delta */
	delta_parameters compute_delta_params(const duration& delta) const;

	/** Convert a matrix to an affine map */
	static math::affine_map<scalar_type> to_map(const matrix_type& M);

	/** Convert an evolution for a given direction to a set of linear constraints
	 *
	 * The evolution is given as a piecewise linear function that is assumed to
	 * be concave.
	 */
	lin_constraint_system<scalar_type> plf_to_constraints(const direction& d,
			const plf_type& f) const;

	/** A functor for measuring the approximation error given a cutpoint */
	std::vector<std::pair<evolution_cache::const_iterator,plif::interval> > evo_with_neighb_intv;


public:
	/** Compute a set of cut points such that all evolutions have a piecewise concave approximation within err distance
	 * of the lower bound such that all inflection points are amidst the cut points.
	 * */
	time_point_set compute_cut_points(
			const cut_point_method& err) const;

	/** Set the relative and absolute error used in numeric computations.
	 *
	 * This is used in the test for inflection points.
	 */
	static void set_numeric_error(error_type num_err);

	/** Get the states at the end of the time domain
	 */
	support_function_provider::ptr get_last_states() const;

	/** Output to stream
	 *
	 * Outputs all of the computed plifs.
	 * The points are separated by std::endl, the values are separated by spaces:
	 * time delta d itv.lower() itv.upper()
	 * The format of the direction and the interval depends on the corresponding
	 * operator<<.
	 */
	void print(std::ostream& os = std::cout) const;

	/** Returns a plf augemented by the error err */
	static plf_type augment_by_error(const plf_type& f, const error_type& err);

private:
	affine_support_flowpipe_description<scalar_type> my_aff;
	time_interval my_time_domain; // time domain
	evolution_cache my_evolutions; // cache of already computed plifs
	delta_parameters_cache_type my_delta_params_cache; //
	delta_parameters_cache_type::const_iterator my_delta_params_it; // iterator to the current entry
	support_function_provider::const_ptr my_Ub; // needed to store the object U+b
	box_type my_S2f;
	cut_point_method my_cut_point_method; // the way to obtain cutpoints
	variable my_time_variable;
	scalar_type min_infeasible_delta; // don't bother computing delta_params higher than this
	scalar_type min_delta; // remember the smallest delta
	bool my_simplify_concave;
	bool my_simplify_convex;

	static error_type my_numeric_error; // = error_type(1e-08,1e-12);
};

}

#include "spacetime_plif.hpp"
#include "spacetime_plif_helpers.hpp"
#include "spacetime_plif_delta_params.hpp"
#include "spacetime_plif_omega.hpp"
#include "spacetime_plif_psi.hpp"
#include "spacetime_plif_convexification.hpp"

#endif /* SPACETIME_PLIF_H_ */
