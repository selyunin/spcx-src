/*
 * finite_hyperbox_utility.h
 *
 *  Created on: Jan 18, 2011
 *      Author: frehse
 */

#ifndef FINITE_HYPERBOX_UTILITY_H_
#define FINITE_HYPERBOX_UTILITY_H_

#include "finite_hyperbox.h"
#include "hyperbox.h"

namespace continuous {

/** Construct the symmetrix bounding box of a support function provider. */
template<typename scalar_type> finite_hyperbox<scalar_type>
finite_bounding_box(const support_function_provider& s);

/** Construct the symmetrix bounding box of a support function provider transformed by an affine map.
 *
 * M is an affine map. */
template<typename scalar_type> finite_hyperbox<scalar_type>
finite_bounding_box(const support_function_provider& s, const math::affine_map<
		scalar_type>& t);

/** Construct the symmetrix bounding box of a finite_hyperbox transformed by an affine map.
 *
 * M is an affine map. */
template<typename scalar_type> finite_hyperbox<scalar_type>
finite_symmetric_bounding_box(const support_function_provider& s,
		const math::affine_map<scalar_type>& t);


/** Construct a general hyperbox from a finite hyperbox */
template<typename scalar_type> hyperbox<scalar_type> construct_hyperbox(
		const finite_hyperbox<scalar_type>& box);

}

#include "finite_hyperbox_utility.hpp"

#endif /* FINITE_HYPERBOX_UTILITY_H_ */
