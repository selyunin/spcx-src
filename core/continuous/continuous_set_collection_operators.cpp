/*
 * continuous_set_collection_operators.cpp
 *
 *  Created on: Mar 30, 2010
 *      Author: frehse
 */

#include "continuous_set_collection_operators.h"
#include "continuous_set_operator_implementations/compute_transformation.h"
#include "continuous_set_transforms/continuous_set_transforms.h"

namespace continuous {

continuous_set::ptr transformer<continuous_set_collection>::compute(
		const continuous_set_collection& c,
		const constant_bound_time_elapse_transform& t) {
	continuous_set_collection::ptr res = continuous_set_collection::ptr(
			new continuous_set_collection());
	for (continuous_set_collection::const_iterator it = c.begin(); it
			!= c.end(); ++it) {
		continuous_set::ptr p = compute_transformation<
				constant_bound_time_elapse_transform> (**it, t);
		res->insert(p);
	}
	return res;
}

continuous_set::ptr transformer<continuous_set_collection>::compute_or_assign(
		continuous_set_collection& c,
		const constant_bound_time_elapse_transform& t) {
	/** @todo one could do an assign to all elements, but
	 * then how to efficiently implement the redundancy checks? */
	/* for (const_iterator it=begin(); it != end(); ++it) {
	 it->swap(it->compute_or_assign_transformation(t));
	 }
	 return get_ptr(); */
	continuous_set_collection::ptr res = continuous_set_collection::ptr(
			new continuous_set_collection());
	for (continuous_set_collection::iterator it = c.begin(); it != c.end(); ++it) {
		res->insert(compute_or_assign_transformation<
				constant_bound_time_elapse_transform> (**it, t));
	}
	return res;
}

continuous_set::ptr transformer<continuous_set_collection>::compute(
		const continuous_set_collection& c, const reset_affine_transform<
				global_types::rational_type>& t) {
	continuous_set_collection::ptr res = continuous_set_collection::ptr(
			new continuous_set_collection());
	for (continuous_set_collection::const_iterator it = c.begin(); it
			!= c.end(); ++it) {
		res->insert(compute_transformation(**it, t));
	}
	return res;
}

continuous_set::ptr transformer<continuous_set_collection>::compute(
		const continuous_set_collection& c, const reset_affine_transform<
				global_types::float_type>& t) {
	continuous_set_collection::ptr res = continuous_set_collection::ptr(
			new continuous_set_collection());
	for (continuous_set_collection::const_iterator it = c.begin(); it
			!= c.end(); ++it) {
		res->insert(compute_transformation(**it, t));
	}
	return res;
}

continuous_set::ptr transformer<continuous_set_collection>::compute_or_assign(
		continuous_set_collection& c, const reset_affine_transform<
				global_types::rational_type>& t) {
	/** @todo one could do an assign to all elements, but
	 * then how to efficiently implement the redundancy checks? */
	/* for (const_iterator it=begin(); it != end(); ++it) {
	 it->swap(it->compute_or_assign_transformation(t));
	 }
	 return get_ptr(); */
	continuous_set_collection::ptr res = continuous_set_collection::ptr(
			new continuous_set_collection());
	for (continuous_set_collection::iterator it = c.begin(); it != c.end(); ++it) {
		res->insert(compute_or_assign_transformation(**it, t));
	}
	return res;
}

continuous_set::ptr transformer<continuous_set_collection>::compute_or_assign(
		continuous_set_collection& c, const reset_affine_transform<
				global_types::float_type>& t) {
	/** @todo one could do an assign to all elements, but
	 * then how to efficiently implement the redundancy checks? */
	/* for (const_iterator it=begin(); it != end(); ++it) {
	 it->swap(it->compute_or_assign_transformation(t));
	 }
	 return get_ptr(); */
	continuous_set_collection::ptr res = continuous_set_collection::ptr(
			new continuous_set_collection());
	for (continuous_set_collection::iterator it = c.begin(); it != c.end(); ++it) {
		res->insert(compute_or_assign_transformation(**it, t));
	}
	return res;
}

continuous_set::ptr transformer<continuous_set_collection>::compute(
		const continuous_set_collection& c, const reset_function_transform& t) {
	continuous_set_collection::ptr res = continuous_set_collection::ptr(
			new continuous_set_collection());
	for (continuous_set_collection::const_iterator it = c.begin(); it
			!= c.end(); ++it) {
		res->insert(compute_transformation(**it, t));
	}
	return res;
}

continuous_set::ptr transformer<continuous_set_collection>::compute_or_assign(
		continuous_set_collection& c, const reset_function_transform& t) {
	/** @todo one could do an assign to all elements, but
	 * then how to efficiently implement the redundancy checks? */
	/* for (const_iterator it=begin(); it != end(); ++it) {
	 it->swap(it->compute_or_assign_transformation(t));
	 }
	 return get_ptr(); */
	continuous_set_collection::ptr res = continuous_set_collection::ptr(
			new continuous_set_collection());
	for (continuous_set_collection::iterator it = c.begin(); it != c.end(); ++it) {
		res->insert(compute_or_assign_transformation(**it, t));
	}
	return res;
}

}
