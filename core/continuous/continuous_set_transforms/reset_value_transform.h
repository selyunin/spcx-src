#ifndef RESET_VALUE_TRANSFORM_H_
#define RESET_VALUE_TRANSFORM_H_

#include "core/continuous/continuous_set_transforms/reset_transform.h"

namespace continuous {

/** Represents the assignment variable := value; where value is of type \p T. */
class reset_value_transform : public reset_transform {
public:
	reset_value_transform(variable_id var, boost::any val) :
		reset_transform(var), my_value(val) {
	}
	;
	virtual ~reset_value_transform() {
	}
	;

	template<typename T> T get_value() {
		return boost::any_cast<T>(my_value);
	}
	;

	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) const {
		used_vars=variable_id_set();
		modif_vars=variable_id_set();
		modif_vars.insert(get_variable_id());
	}
	;

	virtual void print(std::ostream& os) const {
		os << variable(get_variable_id()) << " = " << boost::any_cast<double>(my_value);
	}
	;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;

private:
	boost::any my_value;
};

}

#endif /*RESET_VALUE_TRANSFORM_H_*/
