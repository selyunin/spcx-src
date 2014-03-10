#ifndef CONSTANT_BOUND_TIME_ELAPSE_TRANSFORM_H_
#define CONSTANT_BOUND_TIME_ELAPSE_TRANSFORM_H_

#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"

namespace continuous {

/** Time elapse with derivatives given by constant (state-independent) set p.  */
class constant_bound_time_elapse_transform : public continuous_set_transform {
public:
	typedef boost::shared_ptr<constant_bound_time_elapse_transform> ptr;
	typedef boost::shared_ptr<const constant_bound_time_elapse_transform> const_ptr;

	constant_bound_time_elapse_transform(const continuous_set_const_ptr& p);

	virtual ~constant_bound_time_elapse_transform();

	virtual continuous_set_const_ptr get_relation(continuous_set_const_ptr cset) const;

	virtual continuous_set_const_ptr get_set() const;

	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) const;

	virtual void print(std::ostream& os) const;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;

private:
	continuous_set_const_ptr my_set;
};

}

#endif /*CONSTANT_BOUND_TIME_ELAPSE_TRANSFORM_H_*/
