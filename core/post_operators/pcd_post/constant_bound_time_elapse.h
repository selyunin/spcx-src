#ifndef CONSTANT_BOUND_TIME_ELAPSE_H_
#define CONSTANT_BOUND_TIME_ELAPSE_H_

#include "core/post_operators/continuous_post.h"

namespace hybrid_automata {

/** \brief A continuous post operator for constant_bound_dynamics.
 *
 * It is applicable under the following conditions:
 * - the dynamics are of type constant_bound_dynamics,
 * - the continuous_set class accepts constant_bound_time_elapse_transform.
 */

class constant_bound_time_elapse_post : public continuous_post {
public:
	typedef boost::shared_ptr<constant_bound_time_elapse_post> ptr;
	typedef boost::shared_ptr<const constant_bound_time_elapse_post> const_ptr;

	virtual ~constant_bound_time_elapse_post() {
	}
	;
	virtual constant_bound_time_elapse_post* clone() {
		return new constant_bound_time_elapse_post(*this);
	}
	;

	/** Return the continuous_set that results from applying time elapse
	 * with the time constraints tcons to the continuous set *cset. */
	virtual continuous::continuous_set_ptr post(const time_constraints& tcons,
			const continuous::continuous_set_const_ptr& cset) const;

	static unsigned int max_refine_count;
};

}

#endif /*CONSTANT_BOUND_TIME_ELAPSE_H_*/
