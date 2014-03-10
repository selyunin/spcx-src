#ifndef RESET_TRANSFORM_H_
#define RESET_TRANSFORM_H_

#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"

namespace continuous {

/** Represents the assignment variable := something. Abstract class;
 * derived classes specify what "something" is. */
class reset_transform : public continuous_set_transform {
public:
	reset_transform(variable_id var) :
		my_var_id(var) {
	}
	;
	virtual ~reset_transform() {
	}
	;

	virtual continuous_set_const_ptr get_relation(continuous_set_const_ptr cset) const {
		throw std::runtime_error("no get_relation for reset_transform");
		return continuous_set_const_ptr();
	}
	;
	virtual variable_id get_variable_id() const {
		return my_var_id;
	}
	;

private:
	variable_id my_var_id;
};

}

#endif /*RESET_TRANSFORM_H_*/
