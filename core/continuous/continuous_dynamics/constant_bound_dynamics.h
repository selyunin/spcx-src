#ifndef CONSTANT_BOUND_DYNAMICS_H_
#define CONSTANT_BOUND_DYNAMICS_H_

#include "core/continuous/continuous_dynamics/continuous_dynamics_base.h"

namespace continuous {

/** Time elapse with derivatives given by constant (state-independent) set p.
 *
 * The set of derivatives is stored as a continuous set over the unprimed
 * variables. This is tailored to the time elapse operator of the PPL
 * and avoids unnecessary remapping. */
class constant_bound_dynamics : public continuous_dynamics {
public:
	constant_bound_dynamics(const continuous_set_ptr& p);
	virtual ~constant_bound_dynamics();

	/** Returns the dynamics in predicate form.*/
	virtual dynamics_predicate::ptr get_predicate() const;

	const continuous_set::ptr& get_set() const;

	variable_id_set get_variable_ids() const;

	/** Returns the empty set, to avoid costly projection and difference operations. */
	virtual variable_id_set get_unconstrained_variable_ids() const;

	virtual void print(std::ostream& os) const;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;

private:
	continuous_set::ptr my_set;
};

}

#endif /*CONSTANT_BOUND_DYNAMICS_H_*/
