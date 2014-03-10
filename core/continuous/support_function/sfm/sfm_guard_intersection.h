/*
 * sfm_guard_intersection.h
 *
 *  Created on: May 19, 2010
 *      Author: ray
 */

#ifndef SFM_GUARD_INTERSECTION_H_
#define SFM_GUARD_INTERSECTION_H_

#include <climits>
#include "math/vdom/vdom_vector.h"
#include "math/vdom/lin_constraint.h"
#include "math/numeric/interval.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/continuous/support_function/sfm/sfm_section.h"
#include "core/continuous/support_function/intersection_by_optimization/lb_search_opt.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/support_function_provider_utility.h"

#include <numeric> // for accumulate

//#define _PLOT_INTERSECTION_SAMPLES

namespace continuous{
	namespace support_function{

/*
 * Computes and returns a list of intervals of the indices of the reach sets omega_i's such that
 * each omega_i for i's in an interval intersect with a given guard G and are continuous in the sense that
 * they overlap.
 */

template<class T>
std::list<math::numeric::interval<unsigned int> > get_intersection_interval(
		typename sfm_cont_set<T>::ptr sfm_set_ptr,
		const math::lin_constraint<T>& G) {
	using namespace math::numeric;

	std::list<interval<unsigned int> > g_intersect_intv;

	// First, we treat the degenerate cases

	// If the guard is empty, then it doesn't intersect.
	// so, return an empty list of intervals
	if (math::definitely(!G.is_satisfiable())) {
		LOGGER_OS(DEBUG7,"get_intersection_interval") << "constraint: " << G << " is unsat." << std::endl;
		return g_intersect_intv;
	}
	// If there's a variable in G that's not in the sfm's domain,
	// it always intersects, unless it's unsat (which we just excluded).
	if (math::definitely(G.is_always_satisfied())
			|| !sfm_set_ptr->domain().contains_variables(
					G.get_normal().domain())) {
		interval<unsigned int> intv;
		intv
				= interval<unsigned int> (
						scalar_with_infinity<unsigned int> (0),
						scalar_with_infinity<unsigned int> (
								sfm_set_ptr->get_size() - 1));

		g_intersect_intv.push_back(intv);
		return g_intersect_intv;
	}

	// Now, we treat the normal case

	T G_inh_term = -G.get_canonic_inh_coeff();
	math::vdom_vector<T> G_n = G.get_normal();
	comparison_operator my_sign = G.get_canonic_sign();

	G_n.reorder(sfm_set_ptr->domain());

	// Find the intervals for the constraint
	//    G_n.x <= G_inh_term.
	// We'll deal with equalities later.

	std::pair<unsigned int, T> pr = sfm_set_ptr->extend_sfm(-G_n.get_vector());

	unsigned int neg_index = pr.first;
	T neg_f = pr.second;
//std::cout << "Extended with " << -G_n.get_vector() << " got factor " << neg_f << std::endl;

	const math::matrix<T>& sfm_mat = sfm_set_ptr->get_sfm();
	unsigned int N = sfm_mat.size2();
	unsigned int lw, up;

	bool intv_flag = false; // true if inside guard (inside interval)

	for (unsigned int j = 0; j < N; j++) {
//std::cout << -sfm_mat(neg_index, j) << " < " << G_inh_term * neg_f << " ? ";
		if (math::maybe(is_LE(-sfm_mat(neg_index, j), G_inh_term * neg_f))) {
			if (!intv_flag) {
				intv_flag = true;
				lw = j;
			}
		} else {
			if (intv_flag) {
				up = j - 1;
				interval<unsigned int> intv(lw, up);
				g_intersect_intv.push_back(intv);
				intv_flag = false;
			}
		}
	}
	if (intv_flag) {
		interval<unsigned int> intv(lw, N - 1);
		g_intersect_intv.push_back(intv);
		intv_flag = false;
	}
	if (my_sign == EQ) {
		//std::runtime_error("get_intersection_interval: to be done");

		// we already found the intervals for
		//    G_n.x <= G_inh_term.
		// so now get the ones for
		//    G_n.x >= G_inh_term.
		math::lin_expression<T> lopp(G_n, -G_inh_term);
		math::lin_constraint<T> Gopp(lopp, GE);
		std::list<interval<unsigned int> > g_intersect_intv2 =
				get_intersection_interval(sfm_set_ptr, Gopp);

//		LOGGER_OS(DEBUG7,"get_intersection_interval") << "equality: " << G << " needs to intersects on intervals: " << g_intersect_intv << " and " << g_intersect_intv2 << std::endl;
		g_intersect_intv = intersect(g_intersect_intv, g_intersect_intv2);
	}
	return g_intersect_intv;
}
;

/**
 * Checks if it is necessary to solve the minimization problem to
 * compute the support function of the intersection of the passed
 * continuous set and a guard.
 *
 * @param cont_set_ptr Support function provider set
 * @param G Guard set as a linear constraint
 * @param l Direction at which support function needs to be computed.
 * @param is_bounded Is set to true if the intersection set is a bounded set, false otherwise.
 * @param is_bounded Is set to true if the intersection set is empty, false otherwise.
 * @param to_activate Is set to true if no optimization condition is satisfied.
 * @param max_val Is set to the support function of the intersection of
 * 		  the continuous set and the guard set G if any of the optimization
 * 		  condition holds. to_activate is set to false in this case.
 */
template<typename scalar_type>
void check_opt_conditions(const support_function_provider::const_ptr S,
		const math::lin_constraint<scalar_type>& G,
		const math::vdom_vector<scalar_type>& l_dom,
		bool& is_empty, bool& is_bounded,
		bool& to_activate, scalar_type& max_value){

	double G_inh_term = G.get_canonic_inh_coeff();
	is_bounded = true;
	is_empty = false;

	positional_vdomain set_dom = positional_vdomain(S->get_variable_ids());

	comparison_operator my_sign = G.get_canonic_sign();

	math::lin_expression<double> le = G.get_canonic_l();

	math::vdom_vector<double> G_n = G.get_normal();
	G_n.reorder(set_dom);

	/* Optimization:
	 *
	 */
	math::vdom_vector<double> sup_vec;
	/*Check if the Set S is contained in G
	 */
	if(contains(G,S)){
		S->compute_support(l_dom,max_value,sup_vec,is_empty,is_bounded);
		if(is_empty)
			throw std::runtime_error("sfm_guard_intersection:check_opt_conditions: compute_support requested from EMPTY  set\n");
		if(!is_bounded)
			throw std::runtime_error("sfm_guard_intersection:check_opt_conditions: compute_support requested from UNBOUNDED set\n");
		to_activate = false;
		return;
	}

	/* Check the condition:
	 * 1: G is a hyperplane
	 * 2: l is a multiple of G_n
	 * then return -G_inh_term as the support function value
	 */
	double n_norm = G_n.infinity_norm();
	double l_norm = l_dom.infinity_norm();
	double f = n_norm/l_norm;

	math::vdom_vector<double> d = f * l_dom;
	math::vdom_vector<double> d_neg = -d;

	if (!math::numeric::is_MEQ(f,0.0)){
		if (math::numeric::is_MEQ(G_n.begin(), G_n.end(), d.begin(),d.end())) {// l and n are in the same direction
			max_value = -G_inh_term / f;
			to_activate = false;
			return;
		}
		else if (math::numeric::is_MEQ(G_n.begin(), G_n.end(), d_neg.begin(),d_neg.end())) { // l and n are in the opposite direction
			if(my_sign == EQ){
				max_value =  G_inh_term / f ; // return - inh_term of the guard
				to_activate = false;
				return;
			}
			else{
				S->compute_support(l_dom,max_value,sup_vec,is_empty,is_bounded);
				if(is_empty)
					throw std::runtime_error("sfm_guard_intersection: compute_support requested from EMPTY set\n");
				if(!is_bounded)
					throw std::runtime_error("sfm_guard_intersection: compute_support requested from UNBOUNDED set\n");
				to_activate = false;
				return;
			}
		}
		else{
			// Continue with the lower bound search
			to_activate = true;
		}
	}
	to_activate = true;
};

template<typename scalar_type>
support_function_provider::ptr shift_set(const support_function_provider::const_ptr S,
		const math::vdom_vector<scalar_type>& shift_vector){

	// Create affine map for the translation
	positional_vdomain set_dom = positional_vdomain(S->get_variable_ids());

	math::affine_map<scalar_type> shift_map(set_dom, shift_vector);

	continuous::support_function_provider::ptr shifted_set_ptr(
			new continuous::support_function::sf_unary<scalar_type>(S,
					shift_map));

	return shifted_set_ptr;

};

/** Get the intervals on which an sfm intersects a polyhedron P. */
template<class T>
std::list<math::numeric::interval<unsigned int> > get_intersection_interval(
		typename sfm_cont_set<T>::ptr sfm_set_ptr,
		const continuous::polyhedron<T>& P) {
	typedef math::numeric::interval<unsigned int> intv_type;
	typedef std::list<intv_type> list_type;

	// Initialize the list with the universe
	list_type list;
	if (sfm_set_ptr->get_size() > 0) {
		intv_type I; // universe
		I.set_lower(0);
		I.set_upper(sfm_set_ptr->get_size() - 1);
		list.push_back(I);
	} else {
		// return the empty list
		return list;
	}

	typename math::lin_constraint_system<T>::const_ptr cons =
			P.get_constraints();
	for (typename math::lin_constraint_system<T>::const_iterator i =
			cons->begin(); i != cons->end() && !list.empty(); ++i) {

		// get the list of intervals for constraint *i
		list_type list2 = get_intersection_interval(sfm_set_ptr, *i);

		LOGGER_OS(DEBUG7,"get_intersection_interval") << "constraint: " << *i << " intersects on intervals: " << list2 << std::endl;
		// note that intersect does not return any empty intervals
		list = intersect(list, list2);
	}
	return list;
}
;

template<typename scalar_type>
math::vdom_vector<scalar_type> get_shift_vector(const math::lin_constraint<scalar_type>& G,
		const positional_vdomain set_dom){

	scalar_type G_inh_term = G.get_canonic_inh_coeff();
	scalar_type trans_vec_norm, G_n_norm = 0;

	math::vdom_vector<scalar_type> G_n = G.get_normal();
	G_n.reorder(set_dom);

	for (typename math::vector<scalar_type>::const_iterator it = G_n.begin(); it != G_n.end(); it++) {
		G_n_norm = G_n_norm + (*it) * (*it);
	}

/*
	std::cout << "Shift: Guard:" << G << std::endl;
	std::cout << "Shift: G_n normed:" << G_n_norm << std::endl;
	std::cout << "Shift: G_inh_term:" << G_inh_term;
*/

	math::vdom_vector<scalar_type> trans_vector;

	// @todo I suspect there is a bug here: why divide twice by G_n_norm, and why not taking the sqrt?
	trans_vec_norm =  G_inh_term / G_n_norm;
	trans_vector = trans_vec_norm / G_n_norm * G_n;
	return trans_vector;

}
template<typename scalar_type>
void adjust_shifting(math::numeric::interval<scalar_type>& min_intv,
		const math::vdom_vector<scalar_type>& l,
		const math::vdom_vector<scalar_type>& shift_vector){
	if(min_intv.upper().is_finite())
		min_intv.set_upper(min_intv.upper().get_val() - scalar_product(l, shift_vector));
	if(min_intv.lower().is_finite())
		min_intv.set_lower(min_intv.lower().get_val() - scalar_product(l, shift_vector));

}
/** Compute the support function of the convex hull of the intersection of the omega's
 *  inside the intv, with the guard .
 *
 * @param sfm_ptr pointer to SFM set
 * @param intv Omega interval inside the SFM which intersects with the guard
 * @param l The direction in which to compute the support
 * @param G The Guard as a linear constraint
 *
 * @return The bound on the support function as an interval.
 */
template<typename scalar_type>
typename math::numeric::interval<scalar_type> guard_intersection_support(
		typename sfm_cont_set<scalar_type>::ptr& my_sfm_ptr,math::numeric::interval<unsigned int> intv,
		const math::vdom_vector<scalar_type>& l,
		const math::lin_constraint<scalar_type>& G,
		const std::string minbrak_type,
		double tolerance, unsigned int& samples){

	if(!intv.is_finite()){
		throw std::runtime_error("sfm_guard_intersection: cannot find intersection with infinite sfm section interval");
	}
	LOGGERSW(DEBUG5,"sfm_lin_cons_intersection","Computing support function of sfm section");
	LOGGER_OS(DEBUG5,"sfm_lin_cons_intersection") << intv << ", guard intersection set, @direction: " << l;

	math::numeric::approx_comparator<scalar_type> my_comp;

	unsigned int intv_low = intv.lower().get_val();
	unsigned int intv_up = intv.upper().get_val();

 	std::size_t intv_size = intv.size() + 1;
	std::list<scalar_type> lambdas;  // interval size is guaranteed to be finite here because of previous checks.
	bool is_bounded, is_empty;
	math::vdom_vector<scalar_type> sup_vec;
	scalar_type ext_sup_val;

	std::vector<typename lb_search_opt<scalar_type, sup_functor>::ptr > opt_prbs(intv_size);
	std::vector<bool> relevant_prbs(intv_size,true);
	std::vector<typename support_function_provider::ptr > sfm_section_ptrs(intv_size);

	// everything will be computed on the sfm, so let's take its domain
	positional_vdomain set_dom = my_sfm_ptr->domain(); // positional_vdomain(my_sfm_ptr->get_variable_ids());

	math::vdom_vector<scalar_type> v(l);
	v.reorder(set_dom);

	math::vdom_vector<scalar_type> g_normal(G.get_normal());
	g_normal.reorder(set_dom);

	math::vdom_vector<scalar_type> shift_vector;

	// Obtain the vector to shift the continuous set

	shift_vector = get_shift_vector(G,set_dom);

	scalar_type max, sup_val;
	bool max_set = false, to_activate = true;

	//----
	continuous::support_function::postc_params<scalar_type> post_prb = my_sfm_ptr->get_postc_prb();
	math::affine_map<scalar_type> dynamics_map(set_dom,post_prb.dynamics_A,post_prb.dynamics_b);
	support_function_provider::const_ptr U_ptr = boost::static_pointer_cast<const support_function_provider>(post_prb.input_set_ptr);
	assert(!U_ptr); // Input set has to be a support function provider set.

	//----
	for(int i=intv_low, j=0;i<=intv_up;i++,j++){
		//create omega_i sf provider set as sfm_section(i,i)

		math::numeric::interval<unsigned int> my_intv = math::numeric::interval<unsigned int>(scalar_with_infinity<unsigned int>(i));

		typename support_function_provider::ptr my_sfm_section_ptr(new sfm_section<scalar_type>(my_sfm_ptr, my_intv, false));

		my_sfm_section_ptr = shift_set<scalar_type>(my_sfm_section_ptr, shift_vector);

		//debug: Lets print the set after applying shift

		// Check here if the minimization routine is worth calling for this sfm section member
		check_opt_conditions(my_sfm_section_ptr, G, v, is_empty, is_bounded, to_activate, sup_val);

		if(to_activate) {
			// shift the continuous set.

			// create the functor for the optimization problem

			sup_functor<scalar_type> my_functor(v, g_normal,G.get_inh_coeff(),my_sfm_section_ptr);

			opt_prbs[j] = typename lb_search_opt<scalar_type,sup_functor>::ptr(
					new lb_search_opt<scalar_type, sup_functor>(
							my_functor, minbrak_type, dynamics_map, U_ptr, tolerance));
			sfm_section_ptrs[j] = my_sfm_section_ptr;
		}
		else{
			if(!max_set){
				max = sup_val;
				max_set = true;
			}
			else{
				if(sup_val > max)
					max = sup_val;
			}
			relevant_prbs[j] = false;
		}
		//create optimization prb i;
	}

	const typename sfm_cont_set<scalar_type>::matrix_type& sfm = my_sfm_ptr->get_sfm();
	math::numeric::interval<scalar_type> local_intv;

	bool exit = false;

	math::numeric::interval<scalar_type> omega_interval;
	math::numeric::interval<scalar_type> last_min_interval, min_interval;
	if(max_set){
		last_min_interval.set_lower(max);
		min_interval.set_lower(max);
		last_min_interval.set_upper(max);
		min_interval.set_upper(max);
	}

	while(!exit){
		bool upper_bound_initialised = false;
	    bool lower_bound_initialised = false;
	    //if(min_interval.is_finite())
	    //std::cout << "interval difference:" << min_interval.upper().get_val() - min_interval.lower().get_val() << std::endl;
	    //std::cout << "interval difference:" << min_interval << std::endl;

		for(unsigned int j=0;j<intv_size;j++){
			if(relevant_prbs[j] && opt_prbs[j]){
				// check first if the current bounds on the prb is irrelevant to max interval.
				math::numeric::interval<scalar_type> prb_bounds = opt_prbs[j]->get_min_bounds();
				adjust_shifting(prb_bounds, v, shift_vector);
//				std::cout << "PROBLEM No:" << j << ", min_bounds" << prb_bounds << std::endl;

				if(prb_bounds.upper() < last_min_interval.lower()){
					relevant_prbs[j] = false;
					continue;
				}

				if(prb_bounds.upper().is_finite()){
					if(!upper_bound_initialised){
						upper_bound_initialised = true;
						min_interval.set_upper(prb_bounds.upper().get_val());
					}
					else{
						if(my_comp.is_definitely_strictly_larger(prb_bounds.upper().get_val(), min_interval.upper().get_val()) )
							min_interval.set_upper(prb_bounds.upper().get_val());
					}
				}
				if(prb_bounds.lower().is_finite()){
					if(!lower_bound_initialised){
						lower_bound_initialised = true;
						min_interval.set_lower(prb_bounds.lower().get_val());
					}
					else{
						if(my_comp.is_definitely_strictly_larger(prb_bounds.lower().get_val(), min_interval.lower().get_val()))
							min_interval.set_lower(prb_bounds.lower().get_val());
					}
				}

				if(opt_prbs[j]->get_problem_state() == lb_search_opt<scalar_type, sup_functor>::stop){
//					std::cout << " prob " << j << " reached stop state" <<std::endl;
//					std::cout << "min bounds:" << opt_prbs[j]->get_min_bounds() << std::endl;
					relevant_prbs[j] = false;
					continue;
				}
//				std::cout << "prob no:" << j << ",state:" << opt_prbs[j]->get_problem_state() << std::endl;
				lambdas.push_back(opt_prbs[j]->next_sample());
//				std::cout << "sample request by prob no: " << j <<"=" << opt_prbs[j]->next_sample() << std::endl;
			}
		}

//		std::cout << "lambdas size:" << lambdas.size() << std::endl;
//		std::cout << "min interval is:" << min_interval << std::endl;

		if(lambdas.empty()){
			// all the problem states have reached stop state.
			//checking samples
			unsigned int total_samples = 0;
			for(unsigned int i=intv_low, j=0;i<=intv_up;i++,j++)
				if(opt_prbs[j]){
#ifdef _PLOT_INTERSECTION_SAMPLES
					opt_prbs[j]->plot_graph();
#endif
					total_samples += opt_prbs[j]->get_size();
				}
			samples = total_samples;
			return min_interval;
		}

		typename lb_search_opt<scalar_type,sup_functor>::sample_type new_sample;
		math::vector<scalar_type> ext_dir;
		std::pair<unsigned int, scalar_type> my_pair;
		// update the interval for each of the opt probs and
		// discard the irrelevant optz problems
		//min_interval = math::numeric::interval<scalar_type>();
		for(unsigned int j=intv_low, i=0;j<=intv_up;j++,i++){
			if(relevant_prbs[i] && opt_prbs[i]){

				// we need to make sure that the choosen lambda is within the bounds of the problem
				bool lambda_check = false;
				for(typename std::list<scalar_type>::const_iterator it = lambdas.begin();it!=lambdas.end();it++ ){

					if(!opt_prbs[i]->sample_redundant(*it)){
							new_sample.lambda = *it;
							lambda_check = true;
							break;
					}
				}
				if(!lambda_check){
					relevant_prbs[i] = false;
//					std::cout << "all lambdas redundant for problem no:" << i << std::endl;
					continue;
				}
				// check if all the prbs have become irrelevant here
				bool prb_rel = false;

				for(unsigned int k =0;k<relevant_prbs.size();k++){
					if(relevant_prbs[k]){
						prb_rel = true;
						break;
					}
				}

				if(!prb_rel){

					unsigned int total_samples = 0;
					for(unsigned int i=intv_low, j=0;i<=intv_up;i++,j++)
						if(opt_prbs[j]){
#ifdef _PLOT_INTERSECTION_SAMPLES
								opt_prbs[j]->plot_graph();
#endif
								total_samples += opt_prbs[j]->get_size();
						}
					samples = total_samples;
					return min_interval;
				}

//				std::cout << "choosen lambda:" << new_sample.lambda << "for problem: " << i << std::endl;
				// We need to reorder G's normal before passing it to lb_search_opt.

				ext_dir = lb_search_opt<scalar_type, sup_functor>::map_to_direction(v, g_normal, new_sample.lambda);
				math::vdom_vector<double> l_dom(set_dom, ext_dir);
				sfm_section_ptrs[i]->compute_support(l_dom,ext_sup_val, sup_vec, is_empty,is_bounded);
				new_sample.f_lambda = ext_sup_val;


				opt_prbs[i]->add_sample(new_sample.lambda, new_sample.f_lambda);

				omega_interval = opt_prbs[i]->update_bounds(new_sample);
				if(omega_interval.is_finite() && (my_comp.is_definitely_strictly_smaller(omega_interval.upper().get_val() - omega_interval.lower().get_val(), tolerance)))
					opt_prbs[i]->set_problem_state(lb_search_opt<scalar_type, sup_functor>::stop);

			}
		}

		if(min_interval.is_finite()){
			if(min_interval.upper().get_val() - min_interval.lower().get_val() <= tolerance)
				exit = true;
		}
		lambdas.clear();
	/*	if(update_flag)
			old_interval = min_interval;
	*/
		last_min_interval = min_interval;
	}
	//checking samples
	unsigned int total_samples = 0;
	for(unsigned int i=intv_low, j=0;i<=intv_up;i++,j++)
		if(opt_prbs[j]){
			total_samples += opt_prbs[j]->get_size();
#ifdef _PLOT_INTERSECTION_SAMPLES
			opt_prbs[j]->plot_graph();
#endif
		}
	samples = total_samples;
	//---
	return min_interval;
}
;

/** Computes the median of a list of values
 *
 * Throws if the list is empty.
 */
template<typename scalar_type>
scalar_type compute_median(const std::list<scalar_type>& list) {
	if (list.empty())
		throw std::runtime_error("empty list has no median!");
	// chose first: return *list.begin();

	std::vector<scalar_type> v(list.size());
	std::copy(list.begin(), list.end(), v.begin());
	typename std::vector<scalar_type>::size_type n = v.size() / 2;
	std::nth_element(v.begin(), v.begin() + n, v.end());
	return v[n];
}

/** Computes the average of a list of values
 *
 * Throws if the list is empty.
 */
template<typename scalar_type>
scalar_type compute_average(const std::list<scalar_type>& list) {
	if (list.empty())
		throw std::runtime_error("empty list has no average!");
	// chose first: return *list.begin();

	return std::accumulate(list.begin(), list.end(), scalar_type(0)) / scalar_type(list.size());
}


/**
 * This method solves the lower bound search problem for each of the sfm section member object
 * with the guard constraint simultaneously for optimization purpose, to find the support
 * function of its intersection in the passed direction.
 *
 * @param my_sfm_ptr
 * @param intv
 * @param l
 * @param G
 * @param minbrak_type
 * @param tolerance
 * @param samples
 * @return vector of support function values for each sfm section member
 */
template<typename scalar_type>
typename std::vector<scalar_type> guard_intersection_simult(
		typename sfm_cont_set<scalar_type>::ptr& my_sfm_ptr,math::numeric::interval<unsigned int> intv,
		const math::vdom_vector<scalar_type>& l,
		const math::lin_constraint<scalar_type>& G,
		const std::string minbrak_type,
		double tolerance, unsigned int& samples){

	if(!intv.is_finite()){
		throw std::runtime_error("sfm_guard_intersection_simult: cannot find intersection with infinite sfm section interval");
	}
	LOGGERSW(DEBUG4,"guard_intersection_simult","simultaneous solve");
	LOGGER_OS(DEBUG7,"guard_intersection_simult") << "interval " << intv << ", guard intersection set, @direction:" << l;

	math::numeric::approx_comparator<scalar_type> my_comp;

	unsigned int intv_low = intv.lower().get_val();
	unsigned int intv_up = intv.upper().get_val();

 	std::size_t intv_size = intv.size() + 1;
	std::list<scalar_type> lambdas;  // interval size is guaranteed to be finite here because of previous checks.
	bool is_bounded, is_empty;
	math::vdom_vector<scalar_type> sup_vec;
	scalar_type ext_sup_val;

	std::vector<typename lb_search_opt<scalar_type, sup_functor>::ptr > opt_prbs(intv_size);
	std::vector<bool> active_prbs(intv_size,true);
	std::vector<scalar_type> sf_vector(intv_size); // Contains the result for each opt prb solution.
	std::vector<typename support_function_provider::ptr > sfm_section_ptrs(intv_size);

	positional_vdomain set_dom = my_sfm_ptr->domain();

	math::vdom_vector<scalar_type> v(l);
	v.reorder(set_dom);

	math::vdom_vector<scalar_type> g_normal(G.get_normal());
	g_normal.reorder(set_dom);

	math::vdom_vector<scalar_type> shift_vector;

	// Obtain the vector to shift the continuous set

	//----
	continuous::support_function::postc_params<scalar_type> post_prb = my_sfm_ptr->get_postc_prb();
	math::affine_map<scalar_type> dynamics_map(set_dom,post_prb.dynamics_A,post_prb.dynamics_b);
	support_function_provider::const_ptr U_ptr = boost::static_pointer_cast<const support_function_provider>(post_prb.input_set_ptr);
	assert(!U_ptr); // Input set has to be a support function provider set.
	//----
	shift_vector = get_shift_vector(G,set_dom);

	scalar_type sup_val;
	bool to_activate;

	for(unsigned int i=intv_low, j=0;i<=intv_up;i++,j++){
		//create omega_i sf provider set as sfm_section(i,i)

		math::numeric::interval<unsigned int> my_intv = math::numeric::interval<unsigned int>(scalar_with_infinity<unsigned int>(i));

		typename support_function_provider::ptr my_sfm_section_ptr(new sfm_section<scalar_type>(my_sfm_ptr, my_intv, false));

		my_sfm_section_ptr = shift_set<scalar_type>(my_sfm_section_ptr, shift_vector);

		//debug: Lets print the set after applying shift

		// Check here if the minimization routine is worth calling for this sfm section member
		to_activate = true;
		check_opt_conditions(my_sfm_section_ptr, G, v, is_empty, is_bounded, to_activate, sup_val);

		if(to_activate) {
			// shift the continuous set.
			// create the functor for the optimization problem

			sup_functor<scalar_type> my_functor(v, g_normal,G.get_inh_coeff(),my_sfm_section_ptr);

			opt_prbs[j] = typename lb_search_opt<scalar_type,sup_functor>::ptr(
					new lb_search_opt<scalar_type, sup_functor>(
							my_functor, minbrak_type, dynamics_map, U_ptr, tolerance));
			sfm_section_ptrs[j] = my_sfm_section_ptr;
		}
		else{
			sf_vector[j] = sup_val;
			active_prbs[j] = false;
		}
	}

	// deactivate all prb flags with stop state
	for(unsigned int j=0;j<intv_size;j++){
		if(opt_prbs[j] && opt_prbs[j]->get_problem_state() == lb_search_opt<scalar_type, sup_functor>::stop)
			active_prbs[j] = false;
	}
	bool exit = false;

	math::numeric::interval<scalar_type> omega_interval;

	while(!exit){

	    //1. get the lambda requests for each of the active problems
		scalar_type sample_demand;
		for(unsigned int j=0;j<intv_size;j++){
			if(active_prbs[j] && opt_prbs[j]){

				//std::cout << "prob no:" << j << ",before demand is in state=" << opt_prbs[j]->get_problem_state() << std::endl;
				sample_demand = opt_prbs[j]->next_sample();
				//std::cout << "prob no:" << j << ",after demand is in state=" << opt_prbs[j]->get_problem_state() << std::endl;

				if(opt_prbs[j]->get_problem_state() == lb_search_opt<scalar_type, sup_functor>::stop){
				//		std::cout << " prob " << j << " reached stop state" <<std::endl;
				//		std::cout << "min bounds:" << opt_prbs[j]->get_min_bounds() << std::endl;
						active_prbs[j] = false;
						continue;
				}
				lambdas.push_back(sample_demand);
//				std::cout << "sample request by prob no: " << j <<"=" << sample_demand << std::endl;
			}
		}

//		std::cout << "lambdas size:" << lambdas.size() << std::endl;

		if(lambdas.empty()){
			// all the problem states have reached stop state.
			//checking samples
			unsigned int total_samples = 0;
			for(unsigned int i=intv_low, j=0;i<=intv_up;i++,j++)
				if(opt_prbs[j]){
#ifdef _PLOT_INTERSECTION_SAMPLES
					opt_prbs[j]->plot_graph();
#endif
					total_samples += opt_prbs[j]->get_size();
				}
			samples = total_samples;
			break;
		}

		typename lb_search_opt<scalar_type,sup_functor>::sample_type new_sample;
		math::vector<scalar_type> ext_dir;
		// 2. Choose a lambda from the requested lambdas

		// pick the median of all lambdas as the default value
		//new_sample.lambda = compute_average(lambdas);
		//new_sample.lambda = *lambdas.begin();

		bool lambda_check = false;
		for(unsigned int i=0;i<intv_size;i++){
			if(active_prbs[i] && opt_prbs[i]){

				// we need to make sure that the choosen lambda is within the bounds of the problem
				for(typename std::list<scalar_type>::const_iterator it = lambdas.begin();it!=lambdas.end();it++ ){
					if(!opt_prbs[i]->sample_redundant(*it)){
							new_sample.lambda = *it;
							lambda_check = true;
							break;
					}
				}
				if(!lambda_check){
					std::cout << "guard_intersection_simult: all lambdas redundant for problem no:" << i << std::endl;
					std::cout << "lambdas demand list:" << lambdas << std::endl;
					opt_prbs[i]->print_pivots();
					//opt_prbs[i]->plot_function(-2,1,0.01);
					continue; // Hoping that next iteration will have a non reduntant lambda for this problem
				}

//				// remove lambdas that are redundant for this problem
//				typename std::list<scalar_type>::iterator it = lambdas.begin();
//				while (it != lambdas.end()) {
//					if (opt_prbs[i]->sample_redundant(*it)) {
//						it = lambdas.erase(it);
//					} else {
//						it++;
//					}
//				}
			}
		}
		// GF 2012-12-04 this seems too harsh; there is no reason why all problems have to agree on a single lambda
		if(!lambda_check)
			throw std::runtime_error("all lambdas redundant for all probs! Algorithm flawed!");

//		// if there are nonredundant lambdas left, pick the median of those
//		if(!lambdas.empty())
//			new_sample.lambda = compute_median(lambdas);

//		std::cout << "choosen lambda:" << new_sample.lambda << std::endl;

		// 3. Sample the function for all the problems in the choosen lambda

		ext_dir = lb_search_opt<scalar_type, sup_functor>::map_to_direction(v, g_normal, new_sample.lambda);
		math::vdom_vector<double> l_dom(set_dom, ext_dir);
		for(unsigned int i=0;i<intv_size;i++){

			if(active_prbs[i] && opt_prbs[i]){
//				std::cout << "#prb:" << i << "printing pivots" << std::endl;
//				opt_prbs[i]->print_pivots();


				sfm_section_ptrs[i]->compute_support(l_dom,ext_sup_val, sup_vec, is_empty,is_bounded);
				new_sample.f_lambda = ext_sup_val;
				opt_prbs[i]->add_sample(new_sample.lambda, new_sample.f_lambda);

				omega_interval = opt_prbs[i]->update_bounds(new_sample);
				if(opt_prbs[i]->get_problem_state() == lb_search_opt<scalar_type, sup_functor>::stop){
					active_prbs[i] = false;
				}
				if(omega_interval.is_finite() && (my_comp.is_definitely_strictly_smaller(omega_interval.upper().get_val() - omega_interval.lower().get_val(), tolerance))){
					opt_prbs[i]->set_problem_state(lb_search_opt<scalar_type, sup_functor>::stop);
					active_prbs[i] = false;
				}
			}
		}

		// Check if all the problems are inactive and set exit to true if so.
		exit = true;
		for(unsigned int i=0;i<intv_size;i++){
			if(active_prbs[i]){
				exit = false;
				break;
			}
		}
		//debug
/*		for(unsigned int i=0;i<intv_size;i++){
			if(opt_prbs[i])
				std::cout << "min bounds for #prb:" << i << opt_prbs[i]->get_min_bounds() << std::endl;
		}*/
		//--
		lambdas.clear();
	}
	//checking samples
	unsigned int total_samples = 0;
	for(unsigned int j=0;j<intv_size;j++){
		if(opt_prbs[j]){
			total_samples += opt_prbs[j]->get_size();
#ifdef _PLOT_INTERSECTION_SAMPLES
			opt_prbs[j]->plot_graph();
#endif
		}
		samples = total_samples;
	}
		//---
	// Set the opt prb values with the shift subtracted
	for(unsigned int i=0;i<intv_size;i++){
		if(opt_prbs[i]){
//			std::cout << "min bounds for #prb:" << i << opt_prbs[i]->get_min_bounds() << std::endl;
			sf_vector[i] = opt_prbs[i]->get_min_bounds().upper().get_val() - scalar_product(shift_vector,l);
		}
		else{
			sf_vector[i]-=scalar_product(shift_vector,l);
		}
	}
	return sf_vector;
}
;
/** Compute the support function of the convex hull of the intersection of the omega's
 *  inside the intv, with the guard. This method tries to solve with less samplings.
 *
 * @param sfm_ptr pointer to SFM set
 * @param intv Omega interval inside the SFM which intersects with the guard
 * @param l The direction in which to compute the support
 * @param G The Guard as a linear constraint
 *
 * @return The bound on the support function as an interval.
 */
template<typename scalar_type>
typename math::numeric::interval<scalar_type> guard_intersection_supp_optz(
		typename sfm_cont_set<scalar_type>::ptr& my_sfm_ptr,math::numeric::interval<unsigned int> intv,
		const math::vdom_vector<scalar_type>& l,
		const math::lin_constraint<scalar_type>& G,
		const std::string minbrak_type,
		double tolerance, unsigned int& samples){

	if(!intv.is_finite()){
		throw std::runtime_error("sfm_guard_intersection: cannot find intersection with infinite sfm section interval");
	}
	LOGGERSW(DEBUG5,"guard_intersection_supp_optz","Computing support function of sfm section");
	LOGGER_OS(DEBUG5,"guard_intersection_supp_optz") << intv << ", guard intersection set, @direction: " << l;

	math::numeric::approx_comparator<scalar_type> my_comp;

	unsigned int intv_low = intv.lower().get_val();
	unsigned int intv_up = intv.upper().get_val();

 	std::size_t intv_size = intv.size() + 1;
	std::list<scalar_type> lambdas;  // interval size is guaranteed to be finite here because of previous checks.
	math::vdom_vector<scalar_type> sup_vec;
	scalar_type ext_sup_val;

	std::vector<typename lb_search_opt<scalar_type, sup_functor>::ptr > opt_prbs(intv_size);
	std::vector<bool> relevant_prbs(intv_size,true);
	std::vector<typename support_function_provider::ptr > sfm_section_ptrs(intv_size);

	// since sfm_ptr is const_ptr and we change the sfm object here, lets clone the object
	//typename sfm_cont_set<scalar_type>::ptr my_sfm_ptr = typename sfm_cont_set<scalar_type>::ptr(sfm_ptr->clone());

	positional_vdomain set_dom = my_sfm_ptr->domain();

	math::vdom_vector<scalar_type> v(l);
	v.reorder(set_dom);

	math::vdom_vector<scalar_type> g_normal(G.get_normal());
	g_normal.reorder(set_dom);

	// Obtain the vector to shift the continuous set
	math::vdom_vector<scalar_type> shift_vector;
	shift_vector = get_shift_vector(G,set_dom);

	scalar_type max, sup_val;
	bool is_empty, is_bounded;
	bool max_set = false, to_activate = true;

	continuous::support_function::postc_params<scalar_type> post_prb = my_sfm_ptr->get_postc_prb();
	math::affine_map<scalar_type> dynamics_map(set_dom,post_prb.dynamics_A,post_prb.dynamics_b);
	support_function_provider::const_ptr U_ptr = boost::static_pointer_cast<const support_function_provider>(post_prb.input_set_ptr);
	assert(!U_ptr); // Input set has to be a support function provider set.

	for(unsigned int i=intv_low, j=0;i<= intv_up;i++,j++){
		//create omega_i sf provider set as sfm_section(i,i)

		math::numeric::interval<unsigned int> my_intv = math::numeric::interval<unsigned int>(scalar_with_infinity<unsigned int>(i));

		typename support_function_provider::ptr my_sfm_section_ptr(new sfm_section<scalar_type>(my_sfm_ptr, my_intv, false));

		// Check here if the minimization routine is worth calling for this sfm section member
		check_opt_conditions(my_sfm_section_ptr, G, v, is_empty, is_bounded, to_activate, sup_val);

		if(to_activate) {
			// shift the continuous set.

			my_sfm_section_ptr = shift_set<scalar_type>(my_sfm_section_ptr, shift_vector);
			//debug: Lets print the set after applying shift

			// create the functor for the optimization problem
			sup_functor<scalar_type> my_functor(v, g_normal,G.get_inh_coeff(),my_sfm_section_ptr);

			//debug
/*
			std::ofstream temp_file;
			temp_file.open("/tmp/func_samples.txt");
			for(double i=0; i<=1000;i+=10){
				temp_file << i << " " << my_functor.compute(scalar_type(i)) << std::endl;
			}
			temp_file.close();
			std::system("graph -TX -BC /tmp/func_samples.txt");
*/
			//debug end
			opt_prbs[j] = typename lb_search_opt<scalar_type,sup_functor>::ptr(
					new lb_search_opt<scalar_type, sup_functor>(
							my_functor, minbrak_type, dynamics_map, U_ptr, tolerance));

/*
			// plot the function to minimize
			opt_prbs[j]->plot_function();
*/

			sfm_section_ptrs[j] = my_sfm_section_ptr;
		}
		else{
			if(!max_set){
				max = sup_val;
				max_set = true;
			}
			else{
				if(sup_val > max)
					max = sup_val;
			}
			relevant_prbs[j] = false;
		}
		//create optimization prb i;
	}

	math::numeric::interval<scalar_type> local_intv;

	bool exit = false;

	math::numeric::interval<scalar_type> omega_interval;
	math::numeric::interval<scalar_type> last_min_interval, min_interval;
	if(max_set){
		last_min_interval.set_lower(max);
		min_interval.set_lower(max);
		last_min_interval.set_upper(max);
		min_interval.set_upper(max);
	}

	while(!exit){
		bool upper_bound_initialised = false;
	    bool lower_bound_initialised = false;
	    //if(min_interval.is_finite())
	    //std::cout << "interval difference:" << min_interval.upper().get_val() - min_interval.lower().get_val() << std::endl;
	    //std::cout << "interval difference:" << min_interval << std::endl;

		for(unsigned int j=0;j<intv_size;j++){
			if(relevant_prbs[j] && opt_prbs[j]){
				// check first if the current bounds on the prb is irrelevant to max interval.
				math::numeric::interval<scalar_type> prb_bounds = opt_prbs[j]->get_min_bounds();
				adjust_shifting(prb_bounds, v, shift_vector);
//				std::cout << "PROBLEM No:" << j << ", min_bounds" << prb_bounds << std::endl;

				if(prb_bounds.upper() < last_min_interval.lower()){
					relevant_prbs[j] = false;
					continue;
				}

				if(prb_bounds.upper().is_finite()){
					if(!upper_bound_initialised){
						upper_bound_initialised = true;
						min_interval.set_upper(prb_bounds.upper().get_val());
					}
					else{
						if(my_comp.is_definitely_strictly_larger(prb_bounds.upper().get_val(), min_interval.upper().get_val()) )
							min_interval.set_upper(prb_bounds.upper().get_val());
					}
				}
				if(prb_bounds.lower().is_finite()){
					if(!lower_bound_initialised){
						lower_bound_initialised = true;
						min_interval.set_lower(prb_bounds.lower().get_val());
					}
					else{
						if(my_comp.is_definitely_strictly_larger(prb_bounds.lower().get_val(), min_interval.lower().get_val()))
							min_interval.set_lower(prb_bounds.lower().get_val());
					}
				}

				if(opt_prbs[j]->get_problem_state() == lb_search_opt<scalar_type, sup_functor>::stop){
//					std::cout << " prob " << j << " reached stop state" <<std::endl;
//					std::cout << "min bounds:" << opt_prbs[j]->get_min_bounds() << std::endl;
					relevant_prbs[j] = false;
					continue;
				}
//				std::cout << "prob no:" << j << ",state:" << opt_prbs[j]->get_problem_state() << std::endl;
				lambdas.push_back(opt_prbs[j]->next_sample());
//				std::cout << "sample request by prob no: " << j <<"=" << opt_prbs[j]->next_sample() << std::endl;
			}
		}

//		std::cout << "lambdas size:" << lambdas.size() << std::endl;
//		std::cout << "min interval is:" << min_interval << std::endl;

		if(lambdas.empty()){
			// all the problem states have reached stop state.
			//checking samples
			unsigned int total_samples = 0;
			for(unsigned int i=intv_low, j=0;i<=intv_up;i++,j++)
				if(opt_prbs[j]){
#ifdef _PLOT_INTERSECTION_SAMPLES
					opt_prbs[j]->plot_graph();
#endif
					total_samples += opt_prbs[j]->get_size();
				}
			samples = total_samples;
			return min_interval;
		}

		typename lb_search_opt<scalar_type,sup_functor>::sample_type new_sample;
		math::vector<scalar_type> ext_dir;
		std::pair<unsigned int, scalar_type> my_pair;
		// update the interval for each of the opt probs and
		// discard the irrelevant optz problems
		//min_interval = math::numeric::interval<scalar_type>();
		for(unsigned int j=intv_low, i=0;j<=intv_up;j++,i++){
			if(relevant_prbs[i] && opt_prbs[i]){

				// we need to make sure that the choosen lambda is within the bounds of the problem
				bool lambda_check = false;
				for(typename std::list<scalar_type>::const_iterator it = lambdas.begin();it!=lambdas.end();it++ ){

					if(!opt_prbs[i]->sample_redundant(*it)){
							new_sample.lambda = *it;
							lambda_check = true;
							break;
					}
				}
				if(!lambda_check){
					relevant_prbs[i] = false;
//					std::cout << "all lambdas redundant for problem no:" << i << std::endl;
					continue;
				}
				// check if all the prbs have become irrelevant here
				bool prb_rel = false;

				for(unsigned int k =0;k<relevant_prbs.size();k++){
					if(relevant_prbs[k]){
						prb_rel = true;
						break;
					}
				}

				if(!prb_rel){

					unsigned int total_samples = 0;
					for(unsigned int c=intv_low, a=0;c<=intv_up;c++,a++)
						if(opt_prbs[a]){
#ifdef _PLOT_INTERSECTION_SAMPLES
								opt_prbs[a]->plot_graph();
#endif
								total_samples += opt_prbs[a]->get_size();
						}
					samples = total_samples;
					return min_interval;
				}
//				std::cout << "choosen lambda:" << new_sample.lambda << "for problem: " << i << std::endl;
			}
		}
		ext_dir = lb_search_opt<scalar_type, sup_functor>::map_to_direction(v, g_normal, new_sample.lambda);
		math::vdom_vector<double> l_dom(set_dom, ext_dir);
		for(unsigned int j=intv_low, i=0;j<=intv_up;j++,i++){
			if(opt_prbs[i] && relevant_prbs[i]){
				opt_prbs[i]->add_sample(new_sample.lambda, new_sample.f_lambda);

				sfm_section_ptrs[i]->compute_support(l_dom,ext_sup_val, sup_vec, is_empty,is_bounded);
				new_sample.f_lambda = ext_sup_val;

				omega_interval = opt_prbs[i]->update_bounds(new_sample);
				if(omega_interval.is_finite() && (my_comp.is_definitely_strictly_smaller(omega_interval.upper().get_val() - omega_interval.lower().get_val(), tolerance)))
					opt_prbs[i]->set_problem_state(lb_search_opt<scalar_type, sup_functor>::stop);
			}
		}

		if(min_interval.is_finite()){
			if(min_interval.upper().get_val() - min_interval.lower().get_val() <= tolerance)
				exit = true;
		}
		lambdas.clear();
	/*	if(update_flag)
			old_interval = min_interval;
	*/
		last_min_interval = min_interval;
	} // END of WHILE
	//checking samples
	unsigned int total_samples = 0;
	for(unsigned int i=intv_low, j=0;i<=intv_up;i++,j++)
		if(opt_prbs[j]){
			total_samples += opt_prbs[j]->get_size();
#ifdef _PLOT_INTERSECTION_SAMPLES
			opt_prbs[j]->plot_graph();
#endif
		}
	samples = total_samples;
	//---
	return min_interval;
}
;

} // end of namespace support_function
} // end of namespace continuous

#endif /* SFM_GUARD_INTERSECTION_H_ */
