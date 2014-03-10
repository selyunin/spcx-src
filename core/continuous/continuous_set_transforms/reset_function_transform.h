#ifndef RESET_FUNCTION_TRANSFORM_H_
#define RESET_FUNCTION_TRANSFORM_H_

#include "core/continuous/continuous_set_transforms/reset_transform.h"

namespace continuous {

/** Represents the assignment variable \f$x_k := f(x_1,...,x_n)\f$.
 * For now we use tree_nodes since they are type-indepentent and so can be evaluated
 * with different types during runtime. */
class reset_function_transform : public reset_transform {
public:
	typedef tree::node::ptr function_ptr;
	typedef tree::node::const_ptr function_const_ptr;
	reset_function_transform(variable_id var, const function_ptr& f);
	virtual ~reset_function_transform();

	virtual function_ptr get_function();

	virtual function_const_ptr get_function_const() const;

	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) const;

	virtual void print(std::ostream& os) const;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;

private:
	function_ptr my_function;
};

}

#endif /*RESET_FUNCTION_TRANSFORM_H_*/
