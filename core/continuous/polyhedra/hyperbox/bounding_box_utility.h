/*
 * bounding_box_utility.h
 *
 *  Created on: Apr 6, 2011
 *      Author: frehse
 */

#ifndef BOUNDING_BOX_UTILITY_H_
#define BOUNDING_BOX_UTILITY_H_

#include "bounding_box.h"

namespace continuous {

/** Returns the variables that are zero in all points in a set defined by a support function. */
template<typename scalar_type>
variable_id_set get_vars_bound_to_zero(const support_function_provider& s) {
	hyperbox<scalar_type> bounds = compute_bounding_box<scalar_type> (s);

	// Go through the variables in M_orig and find the ones with bounds zero
	variable_id_set constvars;
	typedef typename hyperbox<scalar_type>::point_type point;
	point lower = bounds.get_l();
	point upper = bounds.get_u();
	const positional_vdomain& dom = bounds.domain();
	for (size_t i = 0; i < dom.size(); ++i) {
		if (lower[i].is_finite() && upper[i].is_finite())
			if (math::numeric::is_MEQ(lower[i].get_val(), scalar_type(0)))
				if (math::numeric::is_MEQ(upper[i].get_val(), scalar_type(0)))
					constvars.insert(dom.get_variable(i).get_id());
	}
	return constvars;
}
;

}

#endif /* BOUNDING_BOX_UTILITY_H_ */
