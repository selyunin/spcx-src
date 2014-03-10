#include "core/continuous/continuous_set_operators.h"

#include <iostream>
#include <stdexcept>

#include "math/tribool.h"

#include "utility/dispatching/double_dispatch.h"

#include "core/continuous/polyhedra/polyhedron_operators.h"

#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"

#include "core/continuous/continuous_set_operator_implementations/intersection.h"
#include "core/continuous/continuous_set_operator_implementations/containment.h"

#include "core/continuous/continuous_set_transforms/sequence_transform.h"
#include "core/continuous/continuous_set_transforms/intersection_transform.h"
#include "core/continuous/continuous_set_transforms/relation_transform.h"

#include "continuous_set_operator_implementations/compute_transformation.h"

#include "core/continuous/support_function/spacetime_flowpipe.h"

namespace continuous {

math::tribool containment_test(const continuous_set* const p1,
		const continuous_set* const p2) {
	// Hack to include containment check for spacetime_flowpipe
	// try casting to flowpipe
	const spacetime_flowpipe<double>* f1 = dynamic_cast<const spacetime_flowpipe<double>*>(p1);
	if (f1) {
		const spacetime_flowpipe<double>* f2 =
				dynamic_cast<const spacetime_flowpipe<double>*>(p2);
		if (f2) {
			return math::definitely(f1->decide_outer_contains(*f2));
		} else {
/*
			std::string name2 = typeid(*p2).name();
			basic_warning("containment_operator",
					"containment_operator not defined for spacetime_flowpipe and "
							+ name2 + ".",
					basic_warning::MISSING_IMPLEMENTATION);
*/
			return math::indeterminate();
		}
	} else {
		const spacetime_flowpipe<double>* f2 =
				dynamic_cast<const spacetime_flowpipe<double>*>(p2);
		if (f2) {
/*			std::string name1 = typeid(*p1).name();
			basic_warning("containment_operator",
					"containment_operator not defined for "
							+ name1 + "and spacetime_flowpipe.",
					basic_warning::MISSING_IMPLEMENTATION);
*/
			return math::indeterminate();
		}
	}

	return dispatching::double_dispatch_tc_cast<math::tribool, containment_operator,
			continuous_set, continuous_set_typelist, containment_test_upcaster>(
			p1, p2);
}


math::tribool containment_test(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2) {
	return containment_test(p1.get(), p2.get());
}

continuous_set_ptr compute_intersection(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2) {
	if (math::definitely(p1->is_universe())) {
		return continuous_set_ptr(p2->clone());
	} else if (math::definitely(p2->is_universe())) {
		return continuous_set_ptr(p1->clone());
	} else {
		return dispatching::double_dispatch_tc<continuous_set_ptr,
				intersection_operator, continuous_set, continuous_set_typelist>(
				p1.get(), p2.get());
	}
}

continuous_set_ptr compute_or_assign_intersection(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2) {
	return compute_intersection(p1, p2);
}

continuous_set_ptr compute_or_assign_union(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2) {
	throw std::runtime_error(
			"missing implementation for compute_or_assign_union");
	return continuous_set_ptr(p1->clone());
}

continuous_set_ptr compute_or_assign_difference(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2) {
	throw std::runtime_error(
			"missing implementation for compute_or_assign_difference");
	return continuous_set_ptr(p1->clone());
}

continuous_set_ptr compute_or_assign_cheap_difference(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2) {
	if (p2->contains(p1))
		return continuous_set_ptr(p1->create_empty());
	else
		return p1;
}

continuous_set_ptr compute_cheap_difference(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2) {
	if (p2->contains(p1))
		return continuous_set_ptr(p1->create_empty());
	else
		return continuous_set_ptr(p1->clone());
}

continuous_set_ptr compute_or_assign_relation_concatenation(
		continuous_set_ptr p1, const continuous_set_const_ptr& p2) {
	// Version without cloning
	p1->reassign_primedness(0, 2); // go from (x,x') to (x'',x')
	p1->reassign_primedness(1, 0); // go from (x'',x') to (x'',x)
	p1 = compute_or_assign_intersection(p1, p2); // results in (x'',x,x')
	p1->existentially_quantify_variables(p1->get_primed_variables(0)); // results in (x'',x')
	p1->reassign_primedness(2, 0); // go from (x'',x') to (x,x')
	return p1;

	/* prettier version but requires cloning
	 continuous_set_ptr q=ps->clone();
	 q->increase_primedness(); // go from (x,x') to (x',x'')
	 intersection_assign(q); // results in (x,x',x'')
	 existentially_quantify_variables(get_primed_variables(1)); // results in (x,x'')
	 decrease_primedness(2); // renames x'' -> x', results in (x,x')
	 */
}

class compute_transformation_visitor: public continuous_set_transform::const_visitor {
public:
	compute_transformation_visitor(const continuous_set& cset) :
		my_res(continuous_set_ptr()), my_cset(cset) {
	}
	;
	continuous_set_ptr get_res() {
		assert(my_res);
		return my_res;
	};
	void dispatch(const intersection_transform* t) {
		my_res = compute_intersection(my_cset.get_const_ptr(), t->get_set());
		assert(my_res);
	}
	;
	void dispatch(const relation_transform* t) {
		my_res = compute_intersection(my_cset.get_const_ptr(), t->get_relation(
				my_cset.get_const_ptr()));
		my_res->existentially_quantify_variables(
				my_res->get_primed_variables(0));
		my_res->decrease_primedness();
		assert(my_res);
	}
	;
	void dispatch(const sequence_transform* t) {
		/* The first transformation is const to create a new continuous_set, the others can
		 * be assigned to this new set. */
		assert(t->begin()!=t->end());
		sequence_transform::const_iterator it = t->begin();
		if (it != t->end()) {
			my_res = compute_transformation(my_cset.get_const_ptr(), *it);
			++it;
		}
		for (; it != t->end(); ++it) {
			my_res = compute_or_assign_transformation(my_res, *it);
		}
		assert(my_res);
	}
	;
	void dispatch(const constant_bound_time_elapse_transform* t) {
		my_res = compute_transformation(my_cset, *t);
		assert(my_res);
	}
	;
	void dispatch(const reset_affine_transform<global_types::rational_type>* t) {
		my_res = compute_transformation(my_cset, *t);
		assert(my_res);
	}
	;
	void dispatch(const reset_affine_transform<global_types::float_type>* t) {
		my_res = compute_transformation(my_cset, *t);
		assert(my_res);
	}
	;
	void dispatch(const reset_function_transform* t) {
		my_res = compute_transformation(my_cset, *t);
		assert(my_res);
	}
	;
	void dispatch(const void* t) {
		const continuous_set_transform* tt=reinterpret_cast<const continuous_set_transform*>(t);
		tt->print(std::cerr);
		std::runtime_error("compute_transformation_visitor: unhandled transform");
	}
	;
private:
	continuous_set::ptr my_res;
	const continuous_set& my_cset;
};

class compute_or_assign_transformation_visitor: public continuous_set_transform::const_visitor {
public:
	compute_or_assign_transformation_visitor(continuous_set& cset) :
		my_res(continuous_set_ptr()), my_cset(cset) {
	}
	;
	continuous_set_ptr get_res() {
		assert(my_res);
		return my_res;
	};
	void dispatch(const intersection_transform* t) {
		my_res
				= compute_or_assign_intersection(my_cset.get_ptr(),
						t->get_set());
		assert(my_res);
	}
	;
	void dispatch(const relation_transform* t) {
		my_res = compute_or_assign_intersection(my_cset.get_ptr(),
				t->get_relation(my_cset.get_const_ptr()));
		my_res->existentially_quantify_variables(
				my_res->get_primed_variables(0));
		my_res->decrease_primedness();
		assert(my_res);
	}
	;
	void dispatch(const sequence_transform* t) {
		assert(t);
		assert(t->begin()!=t->end());
		my_res = my_cset.get_ptr();
		assert(my_res);
		for (sequence_transform::const_iterator it = t->begin(); it != t->end(); ++it) {
			my_res = compute_or_assign_transformation(my_res, *it);
			assert(my_res);
		}
		assert(my_res);
	}
	;
	void dispatch(const constant_bound_time_elapse_transform* t) {
		my_res = compute_or_assign_transformation(my_cset, *t);
		assert(my_res);
	}
	;
	void dispatch(const reset_affine_transform<global_types::rational_type>* t) {
		my_res = compute_or_assign_transformation(my_cset, *t);
		assert(my_res);
	}
	;
	void dispatch(const reset_affine_transform<global_types::float_type>* t) {
		my_res = compute_or_assign_transformation(my_cset, *t);
		assert(my_res);
	}
	;
	void dispatch(const reset_function_transform* t) {
		my_res = compute_or_assign_transformation(my_cset, *t);
		assert(my_res);
	}
	;
	void dispatch(const void* t) {
		const continuous_set_transform* tt=reinterpret_cast<const continuous_set_transform*>(t);
		tt->print(std::cerr);
		std::runtime_error("compute_or_assign_transformation_visitor: unhandled transform");
	}
	;
private:
	continuous_set::ptr my_res;
	continuous_set& my_cset;
};

continuous_set::ptr compute_transformation(continuous_set_const_ptr p,
		const continuous_set_transform_const_ptr& t) {
	assert(t);
	//const_ptr p=get_const_ptr();
	//return t->accept(p);
	compute_transformation_visitor v(*p);
	t->accept(v);
	continuous_set::ptr result=v.get_res();
	assert(result);
	return result;
}

continuous_set::ptr compute_or_assign_transformation(continuous_set_ptr p,
		const continuous_set_transform_const_ptr& t) {
	assert(t);
	//ptr p=get_ptr();
	//return t->accept(p);
	compute_or_assign_transformation_visitor v(*p);
	t->accept(v);
	continuous_set::ptr result=v.get_res();
	assert(result);
	return result;
}

}

