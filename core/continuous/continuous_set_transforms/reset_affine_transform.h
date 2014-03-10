#ifndef RESET_AFFINE_TRANSFORM_H_
#define RESET_AFFINE_TRANSFORM_H_

#include "math/vdom/index_to_variable_id_map_provider.h"
#include "math/vdom/affine_map.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"
#include "core/continuous/predicate_continuous_set.h"

namespace continuous {

/** Represents the affine map \f$x := A * x + b\f$.*/
template<typename scalar_type>
class reset_affine_transform: public continuous_set_transform,
		public math::affine_map<scalar_type> {
public:
	typedef typename math::affine_map<scalar_type>::matrix_type matrix_type;
	typedef typename math::affine_map<scalar_type>::vector_type vector_type;
	using math::affine_map<scalar_type>::compute_used_and_modif_variables;

	/** Convert affine map to transform. */
	reset_affine_transform(const math::affine_map<scalar_type>& M,
			continuous_set::const_ptr U = continuous_set::const_ptr()) :
		math::affine_map<scalar_type>(M), my_U(U) {
		compute_used_and_modif_variables(my_used_vars, my_modif_vars);
	}
	;

	virtual ~reset_affine_transform() {
	}
	;

	virtual continuous_set_const_ptr get_relation(continuous_set_const_ptr cset) const {
		relation_const_ptr p = relation_const_ptr(new predicate_continuous_set(
				this->get_predicate()));
		return p;
	}
	;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const {
		d.dispatch(this);
	}
	;

	/** Obtain the set of used and modified variables */
	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) const {
		used_vars = my_used_vars;
		modif_vars = my_modif_vars;
	}
	;

	/** Obtain the set of nondeterministic inputs */
	virtual const continuous_set::const_ptr& get_input_set() const {
		return my_U;
	}
	;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const {
		math::affine_map<scalar_type>::print(os);
		if (get_input_set()) {
			os << " with offset " << get_input_set();
		}
	}
	;

private:
	variable_id_set my_used_vars;
	variable_id_set my_modif_vars;
	continuous_set::const_ptr my_U;
};

}

#endif /*RESET_AFFINE_TRANSFORM_H_*/
