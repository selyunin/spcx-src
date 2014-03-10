#ifndef CONTINUOUS_SET_TRANSFORM_DECLARATIONS_H_
#define CONTINUOUS_SET_TRANSFORM_DECLARATIONS_H_

#include "utility/typelist.h"
#include "global/global_types.h"

namespace continuous {

/** Register any continuous_set_transform classes in the following typelist. */

/** Forward declaration for inclusion in the typelist. */
class constant_bound_time_elapse_transform;
class intersection_transform;
class relation_transform;
template <typename T> class reset_affine_transform;
class reset_function_transform;
class sequence_transform;

typedef typelist<constant_bound_time_elapse_transform,
typelist<intersection_transform,
typelist<relation_transform,
typelist<reset_affine_transform<global_types::rational_type>,
typelist<reset_affine_transform<global_types::float_type>,
typelist<reset_function_transform,
typelist<sequence_transform,
null_typelist>
> > > > > > continuous_set_transform_typelist;

}

#endif /*CONTINUOUS_SET_TRANSFORM_DECLARATIONS_H_*/
