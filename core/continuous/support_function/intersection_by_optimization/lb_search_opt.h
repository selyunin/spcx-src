/*
 * lb_search_opt.h
 *
 *  Created on: May 7, 2010
 *      Author: ray
 */

#ifndef LB_SEARCH_OPT_H_
#define LB_SEARCH_OPT_H_

#include "convex_optimization.h"
#include "utility/logger_stopwatch.h"
#include <math.h>
#include <iomanip>

namespace continuous{
namespace support_function{

template<class scalar_type, template<typename > class functor>
class lb_search_opt : public convex_opt<scalar_type, functor> {

public:

	typedef global_types::precise_float_type precise_float;
	typedef typename boost::shared_ptr<lb_search_opt<scalar_type, functor> > ptr;
	typedef typename math::numeric::interval<scalar_type> interval;
	typedef typename convex_opt<scalar_type, functor>::min_interval min_interval;
	typedef typename convex_opt<scalar_type, functor>::sample_type sample_type;
	enum STATE{init, one, two, three, four_point, stop};

	lb_search_opt(){};
	lb_search_opt(functor<scalar_type> my_functor, std::string minbrak_type,
			const math::affine_map<scalar_type> dynamics_map,
			const support_function_provider::const_ptr U_ptr,
			const scalar_type eps) : convex_opt<scalar_type, functor>(my_functor, minbrak_type, dynamics_map, U_ptr, eps){
		interval my_interval(scalar_with_infinity<scalar_type>::neg_infty(),
							scalar_with_infinity<scalar_type>::pos_infty());
		math::numeric::approx_comparator<scalar_type> my_comp;
		bool min_found = false;
		if(minbrak_type == "gold_desc")
			//min_found = get_pivots(my_interval,s1,s2,s3,s4);
			minbrak_simple(s1,s2,s3,s4);
		else if(minbrak_type == "parab_desc"){
			minbrak_para_ext(s1,s2,s3,s4);
		}
		else{
			throw std::runtime_error("lb_search_opt: minima bracketing method not known\n");
		}
		if(is_min_found()){
			//std::cout << "min found at min-bracketing\n";
			prob_state = stop;
		}
		else{
			if(my_comp.is_definitely_strictly_larger(s2.f_lambda,s3.f_lambda))
				prob_state = one;
			else if(my_comp.is_definitely_strictly_smaller(s2.f_lambda,s3.f_lambda))
				prob_state = two;
			else prob_state = three;
			min_bounds = interval(scalar_with_infinity<scalar_type>::neg_infty(),
					                    scalar_with_infinity<scalar_type>(get_minimum(s2.f_lambda,s3.f_lambda)));
		}
	};
	lb_search_opt(functor<scalar_type> my_functor,
			const std::string minbrak_type,
			const scalar_type eps) : convex_opt<scalar_type, functor>(my_functor, eps){

		math::numeric::approx_comparator<scalar_type> my_comp;
		bool min_found = false;

		if(minbrak_type == "gold_desc")
			//min_found = get_pivots(intv,s1,s2,s3,s4);
			minbrak_simple(s1,s2,s3,s4);
		else if(minbrak_type == "parab_desc"){
			minbrak_para_ext(s1,s2,s3,s4);
		}
		else{
			throw std::runtime_error("lb_search_opt: minima bracketing method not known\n");
		}

		std::cout << "reached pivots check\n";
		if(is_min_found()){
			prob_state = stop;
			std::cout << "reached pivots check: inside\n";
		}
		else{
			min_bounds = interval(scalar_with_infinity<scalar_type>::neg_infty(),
										scalar_with_infinity<scalar_type>(get_minimum(s2.f_lambda,s3.f_lambda)));

			if(my_comp.is_definitely_strictly_larger(s2.f_lambda,s3.f_lambda))
				prob_state = one;
			else if(my_comp.is_definitely_strictly_smaller(s2.f_lambda,s3.f_lambda))
				prob_state = two;
			else prob_state = three;
		}
	}
	~lb_search_opt(){};

	min_interval  section_search();

	/**
	 * Computes the interval containing the minimum using what we call as lower bound search.
	 */
	min_interval section_search(const interval& search_interval) ;
	/**
	 * Computes the interval containing the minimum using what we call as lower bound search with the initial pivots passed
	 * as arguments.
	 *
	 * @param t1 Pivot 1
	 * @param t2 Pivot 2
	 * @param t3 Pivot 3
	 * @param t4 Pivot 4
	 * @return
	 */
	min_interval section_search(sample_type& t1,sample_type& t2,sample_type& t3,sample_type& t4) ;

	/**
	 * Computes the min containing interval with this classes
	 * new functions, namely update_bounds() and next_sample().
	 * @return minima containing interval.
	 */
	interval section_search_opt();
	/*
	 * Returns the next sampling point request.
	 */
	scalar_type next_sample();

	/**
	 *
	 * @return domain and function values interval containing the function minimum.
	 */
	interval get_min_bounds() const{
		return min_bounds;
	}

	 enum STATE get_problem_state() const {
		 return prob_state;
	 }

	/**
	 * Sets the problem state to the passed state.
	 *
	 * @param s The problem state to be set.
	 */
	 void set_problem_state(enum STATE s) {
		 prob_state = s;
	 }
	/*
	 * Updates the lower and upper bound on the min based on the new sample value.
	 */

	interval update_bounds(sample_type s);

	/*
	 * Returns a lambda value, i.e., sampling point based on the values in the list
	 * passed as argument
	 */
	static scalar_type choose_lambda(std::list<scalar_type>& lambdas);

	/**
	 * Maps parameter lambda to a direction, given by, dir = lambda*c+v
	 *
	 * @param v direction
	 * @param c Normal to the guard
	 * @param lambda parameter
	 * @return direction given by (lambda*c+v)
	 */

	static typename math::vector<scalar_type> map_to_direction(
			const math::vdom_vector<scalar_type>& v,
			const math::vdom_vector<scalar_type>& c,scalar_type lambda);

	bool get_pivots(const interval& my_interval, sample_type& t1,sample_type& t2,sample_type& t3,sample_type& t4);
	/**
	 * Checks if the minima can be deduced from looking at the pivots.
	 * @return true if the minima is found at one of the pivots, false otherwise.
	 */
	bool is_min_found();

	void print_pivots();
	/**
	 *
	 * @return true iff the sample passed is outside the range within which the min is guaranteed to exist.
	 */
	bool sample_redundant(scalar_type lambda);

private:
	enum STATE prob_state;
	sample_type s1, s2, s3, s4, s5; // When the problem is at 4 pivots state, s5 can be made equal to s4.
	interval min_bounds;
};

} // end of namespace support_function
} // end of namespace continuous

#include "lb_search_opt.hpp"

#endif /* LB_SEARCH_OPT_H_ */
