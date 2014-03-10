#ifndef CONTINUOUS_DYNAMICS_DECLARATIONS_H_
#define CONTINUOUS_DYNAMICS_DECLARATIONS_H_

#include "global/global_types.h"
#include "utility/typelist.h"

namespace continuous {

/** Register any continuous_set_transform classes in the following typelist. */

/** Forward declaration for inclusion in the typelist. */
class constant_bound_dynamics;
template <typename T> class ode_affine_dynamics;
class relation_dynamics;

typedef typelist<	constant_bound_dynamics,
		typelist<	ode_affine_dynamics<global_types::rational_type>,
		typelist<	ode_affine_dynamics<global_types::float_type>,
		typelist<	relation_dynamics,
					null_typelist>
		> > > continuous_dynamics_typelist;

}

#endif /*CONTINUOUS_DYNAMICS_DECLARATIONS_H_*/
