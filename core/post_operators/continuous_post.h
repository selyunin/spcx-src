/*
 * continuous_post.h
 *
 *  Created on: Jun 23, 2009
 *      Author: frehse
 */

#ifndef CONTINUOUS_POST_H_
#define CONTINUOUS_POST_H_

#include "core/post_operators/post_operator.h"

/** Forward declaration of classes used in header file. */
namespace discrete {
class discrete_set;
typedef boost::shared_ptr<discrete_set> discrete_set_ptr;
typedef boost::shared_ptr<const discrete_set> discrete_set_const_ptr;
}
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
}
namespace hybrid_automata {
class time_constraints;
}

namespace hybrid_automata {

/** A continuous post operator is a post operator that iterates
 * over a set of discrete sets provided by the automaton.
 *
 * For a full implementation, it suffices to implement clone()
 * and post(time_constraints,continuous_set).
 */

class continuous_post: public post_operator {
public:
	typedef boost::shared_ptr<continuous_post> ptr;
	typedef boost::shared_ptr<const continuous_post> const_ptr;

	/** Default constructor.
	 *
	 *
	 * It's main point is to define parameter values that make sense
	 * (rather than, say, a time horizon of zero).
	 * */
	continuous_post() : my_time_horizon(-1.0), my_sampling_time(-1.0) {};

	virtual ~continuous_post() {
	}
	;
	virtual continuous_post* clone() = 0;

	/** Return the continuous_set that results from applying time elapse
	 * with the time constraints tcons to the continuous set *cset.
	 *
	 * If the returned set is empty, the implementation must return
	 * a null pointer (continuous_set_ptr()). This is to avoid
	 * emptiness checks further down the line.
	 *
	 * This is the only function (apart from clone()) that needs to be provided by derived
	 * implementations.*/
	virtual continuous::continuous_set_ptr post(const time_constraints& tcons,
			const continuous::continuous_set_const_ptr& cset) const = 0;

	/** Return the continuous_set that results from applying time elapse
	 * in the discrete set *dset (all locations of which are considered
	 * equivalent w.r.t. time elapse) to the continuous set *cset.
	 *
	 * The default implementation simply chooses the time constraints of
	 * the first location in the set.*/
	virtual continuous::continuous_set_ptr post(
			const hybrid_automaton_const_ptr& aut,
			const discrete::discrete_set_const_ptr& dset,
			const continuous::continuous_set_const_ptr& cset) const;

	/** Apply the post operator to the symbolic state sstate and
	 * add the result to both result sets. */
	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const symbolic_state_ptr& sstate) const;

	using post_operator::add_post_states;
	/** Set time horizon to th.
	 *
	 * A negative value denotes unlimited time. */
	virtual void set_time_horizon(double th);
	virtual double get_time_horizon() const;
	/** Set sampling time to ts.
	 *
	 * A negative value denotes no sampling. */
	virtual void set_sampling_time(double ts);
	virtual double get_sampling_time() const;

private:
	double my_time_horizon;
	double my_sampling_time;
};

}

#endif /* CONTINUOUS_POST_H_ */
