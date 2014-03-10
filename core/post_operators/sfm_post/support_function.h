#ifndef SUPPORT_FUNCTION_H_
#define SUPPORT_FUNCTION_H_

#include "math/matrix.h"
#include "math/matrix_operators.h"
#include "math/vector.h"
#include "math/scalar_types/scalar_with_infinity.h"
#include "core/continuous/polyhedra/polyhedron.h"
#include "core/continuous/polyhedra/polyhedron_utility.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/ppl_polyhedron/convert_to_ppl.h"
#include "core/discrete/discrete_set_stl_set.h"
#include "core/continuous/polyhedra/hyperbox/bounding_box.h"
#include "core/continuous/support_function/template_directions/choose_directions.h"

#include "omega_model_factory.h"

namespace continuous {
namespace support_function {

/**
 * Computes the set of symbolic states that is reachable from cset with elapse of time upto delta.
 *
 * \param csup The initial set from which the time elapse reachable set is computed
 * \param M The affine dynamics x'=Ax+b.
 * \param inv The invariant associated with the location.
 * \param time_step The time step.
 * \param time_horizon The time bound.
 *
 * Returns the set of symbolic states which contains the reachable time elapse set from cset. By symbolic state,
 * we mean a pair of discrete set and a continuous set. sstate_list is a collection of such symbolic states in which
 * the discrete set contains only one integer element k and the continuous set is a subset of R^n which is reachable
 * in the kth iteration of the discrete time elapse set computation algorithm. Here, the discrete time elapse computation
 * algorithm is the discrete version of the reachability algorithm using support function [Girard et al].
 */
template<class scalar_type> typename sfm_cont_set<scalar_type>::ptr
support_function_time_elapse(const support_function_provider::const_ptr& csup,
		const math::affine_map<scalar_type>& M, support_function_provider::const_ptr U, const typename polyhedron<
				scalar_type>::const_ptr& inv, double time_horizon,
		double time_step, double step_tolerance = -1.0);
}
}

#include "support_function.hpp"
#endif /*SUPPORT_FUNCTION_H_*/
