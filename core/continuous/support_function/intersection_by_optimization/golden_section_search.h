#ifndef _GOLDEN_SECTION_SEARCH
#define _GOLDEN_SECTION_SEARCH

#include "convex_optimization.h"
#include "utility/logger_stopwatch.h"

namespace continuous {
namespace support_function{

/**
 * This function is to provide interfaces required to compute precise intersection operation between
 * a convex set and a hyperplane c.x=d using golden section search algorithm.
 *
 * @author Rajarshi
 */

template<class scalar_type, template<typename > class functor>
class golden_section_search : public convex_opt<scalar_type, functor >{
public:

	typedef typename math::numeric::interval<scalar_type> interval;
	typedef typename convex_opt<scalar_type, functor>::min_interval min_interval;

	/**
	 * This structure is a function and argument pair, where the function from Reals to Reals.
	 */

	/* Constructors */
	golden_section_search(){};

	golden_section_search(functor<scalar_type> my_functor, const scalar_type eps) :
		convex_opt<scalar_type, functor>(my_functor, eps){}
	~golden_section_search(){};

	min_interval section_search();

	min_interval section_search(const interval& search_interval);
};

} // end of namespace support_function
} // end of namespace continuous

#include "golden_section_search.hpp"

#endif /* _GOLDEN_SECTION_SEARCH */


