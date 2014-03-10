/*
 * lower_bound_search.h
 *
 *  Created on: Feb 21, 2010
 *      Author: ray
 */

#ifndef LOWER_BOUND_SEARCH_H_
#define LOWER_BOUND_SEARCH_H_

#include "convex_optimization.h"
#include "utility/logger_stopwatch.h"
#include <math.h>
#include <climits>

//#define OPT_INIT_SAMPLE__

/**
 * Convex Optimization using lower bound search technique
 *
 * @author Rajarshi
 */

/** Switches on plotting of samples taken during the intersection using optimization. */
//#define PLOT_INTERSECTION_SAMPLES

namespace continuous {
namespace support_function {

template<class scalar_type, template<typename > class functor>
class lower_bound_search : public convex_opt<scalar_type, functor> {
public:
	typedef global_types::precise_float_type precise_float;
	typedef typename math::numeric::interval<scalar_type> interval;
	typedef typename convex_opt<scalar_type, functor>::sample_type sample_type;
	typedef typename convex_opt<scalar_type, functor>::min_interval min_interval;
	/* Constructors */
	lower_bound_search(){};

	lower_bound_search(functor<scalar_type> my_functor, const std::string minbrak_type,
			const math::affine_map<scalar_type> dynamics_map,
			const support_function_provider::const_ptr U_ptr,
			const scalar_type eps) : convex_opt<scalar_type, functor>(my_functor, minbrak_type, dynamics_map, U_ptr, eps){}

	lower_bound_search(functor<scalar_type> my_functor, const std::string minbrak_type,
				const scalar_type eps) : convex_opt<scalar_type, functor>(my_functor, minbrak_type , eps){}
	~lower_bound_search(){};

	min_interval  section_search();
	/**
	 * Computes the interval containing the minimum using what we call as lower bound search.
	 */
	min_interval section_search(const interval& search_interval) ;
	/**
	 * Search is bounded on the number of function samples.
	 *
	 * @param search_interval
	 * @param bound The sample bound
	 * @return
	 */
	min_interval section_search_bounded(const interval& search_interval,const unsigned int bound) ;
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
	 * Search bounded on the number of function samples.
	 *
	 * @param t1
	 * @param t2
	 * @param t3getting max double value
	 * @param t4
	 * @param bound The bound on the function samples.
	 * @return
	 */
	min_interval section_search_bounded(sample_type& t1,sample_type& t2,sample_type& t3,sample_type& t4, unsigned int bound) ;
	/**
	 * Computes the 4 initial pivot points to start with. The pivot points should satisfy the
	 * following condition:
	 * s2.f_lambda <= s1.f_lambda and s3.f_lambda <= s4.f_lambda
	 *
	 * @return true if 4 initial pivots are found successfully, false otherwise
	 */
	bool get_pivots(const interval& my_interval, sample_type& t1,sample_type& t2,sample_type& t3,sample_type& t4, const unsigned int bound) ;
	/**
	 * Checks if the minima can be deduced from looking at the pivots.
	 *
	 * @param my_interval
	 * @param s1
	 * @param s2
	 * @param s3
	 * @param s4
	 * @return true if the minima is found at one of the pivots, false otherwise.
	 */
	bool is_min_found(sample_type& s1,sample_type& s2,sample_type& s3,sample_type& s4);
	void print_pivots(const sample_type& s1, const sample_type& s2, const sample_type& s3, const sample_type& s4, const sample_type& s5);

};

} // end of namespace support_function
} // end of namespace continuous

#include "lower_bound_search.hpp"


#endif /* LOWER_BOUND_SEARCH_H_ */
