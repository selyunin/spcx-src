/*
 * direct_discrete_post.h
 *
 *  Created on: Aug 19, 2009
 *      Author: frehse
 */

#ifndef DIRECT_DISCRETE_POST_H_
#define DIRECT_DISCRETE_POST_H_

#include "core/post_operators/discrete_post.h"

namespace hybrid_automata {

/** \brief A discrete post operator for any transform that is accepted by
 * the continuous_set class.
 *
 * It is applicable under the following conditions:
 * - the continuous_set class accepts the intersection
 * - the continuous_set class accepts the transform type of the transition.
 */

class direct_discrete_post: public discrete_post {
public:
	typedef boost::shared_ptr<direct_discrete_post> ptr;
	typedef boost::shared_ptr<const direct_discrete_post> const_ptr;

	virtual ~direct_discrete_post();
	virtual direct_discrete_post* clone();

	/** Return the continuous_set that results from applying the transition
	 * *trans to the continuous set *cset. */
	virtual continuous::continuous_set_collection post(const jump_constraints& trans,
			continuous::continuous_set_const_ptr source_inv,
			continuous::continuous_set_const_ptr target_inv,
			continuous::continuous_set_const_ptr cset) const;

};

}

#endif /* DIRECT_DISCRETE_POST_H_ */
