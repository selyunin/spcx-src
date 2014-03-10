#ifndef RELATION_DYNAMICS_H_
#define RELATION_DYNAMICS_H_

#include "core/continuous/continuous_dynamics/continuous_dynamics_base.h"

namespace continuous {

/** Represents the relation \f$ \{ (x,x') \mid x' \in F(x) \}\f$.
 * The relation is represented as
 * a \p continous_set, in which the unprimed variables represent x and the
 * (singly-) primed variables represent x'.  */
class relation_dynamics : public continuous_dynamics {
public:
	relation_dynamics(const relation::ptr& p);

	virtual ~relation_dynamics();

	virtual dynamics_predicate::ptr get_predicate() const;

	variable_id_set get_variable_ids() const;

	/** Returns the empty set, to avoid costly projection and difference operations. */
	virtual variable_id_set get_unconstrained_variable_ids() const;

	virtual void print(std::ostream& os) const;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;

	virtual relation::const_ptr get_relation() const;
	virtual void set_relation(relation::ptr p);

private:
	relation::ptr my_relation;
};

}

#endif /*RELATION_DYNAMICS_H_*/
