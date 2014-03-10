#ifndef CONTINUOUS_SET_TRANSFORM_COMPOSITION_H_
#define CONTINUOUS_SET_TRANSFORM_COMPOSITION_H_

#include "boost/shared_ptr.hpp"

#include "math/vdom/variable.h"

/** Forward declarations */
namespace continuous {
class continuous_set_transform;
typedef boost::shared_ptr<continuous_set_transform> continuous_set_transform_ptr;
typedef boost::shared_ptr<const continuous_set_transform> continuous_set_transform_const_ptr;
}

namespace continuous {

/** Computes the parallel composition of the two transformations t1 and t2, and puts the result in t_ret.
 * Returns true iff the composition was successful. */
bool
parallel_compose(continuous_set_transform_ptr& t_ret,
		const continuous_set_transform_const_ptr& t1,
		const continuous_set_transform_const_ptr& t2);

continuous_set_transform_ptr parallel_compose(const continuous_set_transform_const_ptr& t1,
		const continuous_set_transform_const_ptr& t2);

/** Extend domain and codomain of the transform such that the variables in vis remain constant.
 * If no such transform exists, return a null pointer.
 */
continuous_set_transform_ptr extend_transform_with_const_variables(
		const continuous_set_transform_const_ptr& t1, const variable_id_set& vis);

}

#endif /*CONTINUOUS_SET_TRANSFORM_COMPOSITION_H_*/
