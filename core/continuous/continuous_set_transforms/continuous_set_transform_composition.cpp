#include "core/continuous/continuous_set_transforms/continuous_set_transform_composition.h"
//#include "../continuous_set.h"

#include "core/continuous/continuous_set_transforms/continuous_set_transforms.h"
#include "core/continuous/continuous_set_operators.h"
//#include "continuous_set_transforms.h"

#include "core/continuous/polyhedra/polyhedron.h"
#include "core/continuous/polyhedra/polyhedron_utility.h"

#include "math/vdom/affine_map_utility.h"
#include "core/predicates/valuation_function_tree_utility.h"

#include "utility/basic_warning.h"

namespace continuous {

/** Declaration of the parallel composition operator and its wrapper. */
template<typename T1, typename T2> class parallel_composition_operator {
public:
	static continuous_set_transform::ptr implement(const T1* t1, const T2* t2) {
		assert(t1);
		assert(t2);
		continuous_set_transform::ptr t_ret;
		// try to sequentialize them
		variable_id_set used1, modif1, used2, modif2;
		t1->get_used_and_modif_variables(used1, modif1);
		t2->get_used_and_modif_variables(used2, modif2);
		if (set_is_disjoint<variable_id> (used1, modif2)) {
			t_ret = continuous_set_transform::ptr(new sequence_transform(t2->get_const_ptr(),
					t1->get_const_ptr()));
			return t_ret;
		} else if (set_is_disjoint(used2, modif1)) {
			t_ret = continuous_set_transform::ptr(new sequence_transform(t1->get_const_ptr(),
					t2->get_const_ptr()));
			return t_ret;
		}
		return continuous_set_transform::ptr(); // failed to find composition
	}
	;
};

bool parallel_compose(continuous_set_transform::ptr& t_ret,
		const continuous_set_transform::const_ptr& t1,
		const continuous_set_transform::const_ptr& t2) {
	assert(t1);
	assert(t2);
	t_ret = parallel_compose(t1, t2);
	return t_ret;
}

continuous_set_transform::ptr parallel_compose(const continuous_set_transform::const_ptr& t1,
		const continuous_set_transform::const_ptr& t2) {
	assert(t1);
	assert(t2);
	return dispatching::double_dispatch_tc<continuous_set_transform::ptr,
			parallel_composition_operator, continuous_set_transform,
			continuous_set_transform_typelist>(t1.get(), t2.get());
}

class extend_transform_with_const_variables_visitor: public continuous_set_transform::const_visitor {
public:
	extend_transform_with_const_variables_visitor(const variable_id_set& vis) :
		my_vis(vis), my_ptr(continuous_set_transform::ptr()) {
	}
	;
	void dispatch(const constant_bound_time_elapse_transform* t) {
		throw std::runtime_error("missing implementation in extend_transform_with_const_variables_visitor");
	}
	;
	void dispatch(const intersection_transform* t) {
		throw std::runtime_error("missing implementation in extend_transform_with_const_variables_visitor");
	}
	;
	void dispatch(const relation_transform* t) {
		continuous_set::ptr cset(t->get_relation_const()->clone());
		if (polyhedron<Rational>::ptr poly=boost::dynamic_pointer_cast<polyhedron<Rational> >(cset)) {
			// Add a const constraint for each variable in vis
			add_const_relation_constraints(*poly, my_vis);
			my_ptr=continuous_set_transform::ptr(new relation_transform(poly));
			//std::cerr << "extended: "<< *poly;
		}
		else if(predicate_continuous_set::ptr pred_cset = boost::dynamic_pointer_cast<predicate_continuous_set> (cset)){
			predicate_continuous_set::predicate_type_ptr pred = pred_cset->get_predicate();
			valuation_functions::add_const_relation_constraints(pred, my_vis);
			pred_cset->set_predicate(pred);
			my_ptr=continuous_set_transform::ptr(new relation_transform(pred_cset));
		}
		else {
			throw std::runtime_error("unknown relation_transform in extend_transform_with_const_variables_visitor");
		}
	}
	;
	void dispatch(const reset_affine_transform<global_types::rational_type>* t) {
		positional_vdomain dom(my_vis);
		math::affine_map<global_types::rational_type> const_map(dom,global_types::rational_type(0));
		try {
			math::info_resolver<global_types::rational_type> resolv;
			math::affine_map<global_types::rational_type> new_map = compose(
					(math::affine_map<global_types::rational_type>) (*t),
					const_map, resolv);
			if (resolv.has_conflict()) {
				basic_warning("Composing transform", "Created empty transform",
						basic_warning::UNUSUAL_INPUT);
				positional_vdomain cdom = compose(t->domain(), dom);
				new_map = math::affine_map<global_types::rational_type>::void_map(cdom);
			}
			my_ptr = continuous_set_transform::ptr(new reset_affine_transform<
					global_types::rational_type>(new_map));
		} catch (std::exception& e) {
			std::stringstream sdyn,svars;
			logger::copyfmt_to(svars);
			sdyn << *t;
			print_variable_id_set(svars,my_vis);
			throw basic_exception("Could not enforce on the following dynamics that variables " +svars.str()+" remain constant:\n"+sdyn.str(),e);
		}
	}
	;
	void dispatch(const reset_affine_transform<global_types::float_type>* t) {
		positional_vdomain dom(my_vis);
		math::affine_map<global_types::float_type> const_map(dom,global_types::float_type(0));
		try {
			math::info_resolver<global_types::float_type> resolv;
			math::affine_map<global_types::float_type> new_map = compose(
					(math::affine_map<global_types::float_type>) (*t),
					const_map, resolv);
			if (resolv.has_conflict()) {
				basic_warning("Composing transform", "Created empty transform",
						basic_warning::UNUSUAL_INPUT);
				new_map = math::affine_map<global_types::float_type>();
			}
			my_ptr = continuous_set_transform::ptr(new reset_affine_transform<
					global_types::float_type>(new_map));
		} catch (std::exception& e) {
			std::stringstream sdyn,svars;
			logger::copyfmt_to(svars);
			sdyn << *t;
			print_variable_id_set(svars,my_vis);
			throw basic_exception("Could not enforce on the following transform that variables " +svars.str()+" remain constant:\n"+sdyn.str(),e);
		}
	}
	;
	void dispatch(const reset_function_transform* t) {
		throw std::runtime_error("missing implementation in extend_transform_with_const_variables_visitor");
	}
	;
	void dispatch(const sequence_transform* t) {
		throw std::runtime_error("missing implementation in extend_transform_with_const_variables_visitor");
	}
	;
	continuous_set_transform::ptr get_transform() const {
		return my_ptr;
	}
	;
private:
	const variable_id_set& my_vis;
	continuous_set_transform::ptr my_ptr;
};

continuous_set_transform::ptr extend_transform_with_const_variables(
		const continuous_set_transform::const_ptr& t1, const variable_id_set& vis) {
	if (t1) {
		extend_transform_with_const_variables_visitor v(vis);
		t1->accept(v);
		return v.get_transform();
	}
	else
		return continuous_set_transform::ptr();
}

//--------------------------------------------------------------------------

/** Compose two reset_affine_transforms */
/** Compose a relation with a relation. Needed to avoid ambiguities. */
template<typename scalar_type> class parallel_composition_operator<reset_affine_transform<scalar_type>, reset_affine_transform<scalar_type> > {
public:
	static continuous_set_transform::ptr implement(const reset_affine_transform<scalar_type>* t1,
			const reset_affine_transform<scalar_type>* t2) {
		assert(t1);
		assert(t2);
		// compose the transforms, throw if not possible
		// use compose operator of affine_map
		math::affine_map<scalar_type> new_map;
		math::info_resolver<scalar_type> resolv;
		new_map = math::compose<scalar_type>(*t1, *t2, resolv);
		if (resolv.has_conflict()) {
			basic_warning("Composing transform", "Created empty transform",
					basic_warning::UNUSUAL_INPUT);
			positional_vdomain cdom = compose(t1->domain(), t2->domain());
			new_map = math::affine_map<scalar_type>::void_map(cdom);
		}
		continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(
				new reset_affine_transform<scalar_type> (new_map));
		return t_ret;
	}
	;
};

continuous_set_transform::ptr get_intersection(const intersection_transform* t1,
		const intersection_transform* t2);

continuous_set_transform::ptr get_intersection(const pre_intersection_transform* t1,
		const pre_intersection_transform* t2);

continuous_set_transform::ptr get_intersection(const post_intersection_transform* t1,
		const post_intersection_transform* t2);

continuous_set_transform::ptr parallel_compose_sequence(const sequence_transform* t1,
		const sequence_transform* t2, bool try_swapping = true);

/** Compose a sequence with a sequence. */
template<> class parallel_composition_operator<sequence_transform, sequence_transform> {
public:
	static continuous_set_transform::ptr implement(const sequence_transform* t1,
			const sequence_transform* t2) {
		assert(t1);
		assert(t2);
		return parallel_compose_sequence(t1, t2);
	}
	;
};

/** Compose a sequence with a transform. */
template<typename T2> class parallel_composition_operator<sequence_transform, T2> {
public:
	static continuous_set_transform::ptr implement(const sequence_transform* t1, const T2* t2) {
		assert(t1);
		assert(t2);
		// create a singleton sequence
		sequence_transform::ptr t3 = sequence_transform::ptr(new sequence_transform(
				t2->get_const_ptr()));
		return parallel_composition_operator<sequence_transform, sequence_transform>::implement(t1,
				t3.get());
	}
	;
};

/** Compose a transform with a sequence. */
template<typename T1> class parallel_composition_operator<T1, sequence_transform> {
public:
	static continuous_set_transform::ptr implement(const T1* t1, const sequence_transform* t2) {
		assert(t1);
		assert(t2);
		return parallel_composition_operator<sequence_transform, T1>::implement(t2, t1);
	}
	;
};

/** Compose a relation with a transform. */
template<typename T> class parallel_composition_operator<relation_transform, T> {
public:
	static continuous_set_transform::ptr implement(const relation_transform* t1, const T* t2) {
		assert(t1);
		assert(t2);
		// intersect the relations
		continuous_set_ptr cset = compute_intersection(t1->get_relation_const(), t2->get_relation(
				t1->get_relation_const()));
		continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(new relation_transform(
				cset));
		return t_ret;
	}
	;
};

/** Compose a transform with a relation. */
template<typename T> class parallel_composition_operator<T, relation_transform> {
public:
	static continuous_set_transform::ptr implement(const T* t1, const relation_transform* t2) {
		assert(t1);
		assert(t2);
		return parallel_composition_operator<relation_transform, T>::implement(t2, t1);
	}
	;
};

/** Compose a relation with a relation. Needed to avoid ambiguities. */
template<> class parallel_composition_operator<relation_transform, relation_transform> {
public:
	static continuous_set_transform::ptr implement(const relation_transform* t1,
			const relation_transform* t2) {
		assert(t1);
		assert(t2);
		// intersect the relations
		continuous_set_ptr cset = continuous_set_ptr();
		cset = compute_intersection(t1->get_relation_const(), t2->get_relation_const());
		continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(new relation_transform(
				cset));

		// if in DEBUG, show result
		IFLOGGER(DEBUG7) {
			LOGGER_OS(DEBUG7,__FUNCTION__) << "composing transition relation " << t1->get_relation_const() << " and " << t2->get_relation_const() << ", result: " << cset;
		}
		return t_ret;
	}
	;
};

/** Compose a relation with a sequence. Needed to avoid ambiguities. */
template<> class parallel_composition_operator<relation_transform, sequence_transform> {
public:
	static continuous_set_transform::ptr implement(const relation_transform* t1,
			const sequence_transform* t2) {
		assert(t1);
		assert(t2);
		// intersect the relations
		continuous_set_ptr cset = continuous_set_ptr();
		cset = compute_intersection(t1->get_relation_const(), t2->get_relation(
				t1->get_relation_const()));
		continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(new relation_transform(
				cset));
		return t_ret;
	}
	;
};

/** Compose a relation with a sequence. Needed to avoid ambiguities. */
template<> class parallel_composition_operator<sequence_transform, relation_transform> {
public:
	static continuous_set_transform::ptr implement(const sequence_transform* t1,
			const relation_transform* t2) {
		assert(t1);
		assert(t2);
		return parallel_composition_operator<relation_transform, sequence_transform>::implement(t2,
				t1);
	}
	;
};

/** Compose a pre_intersection_transform with a transform. */
template<typename T> class parallel_composition_operator<pre_intersection_transform, T> {
public:
	static continuous_set_transform::ptr implement(const pre_intersection_transform* t1,
			const T* t2) {
		assert(t1);
		assert(t2);
		// t1 goes in front
		continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(new sequence_transform(
				t1, t2));
		return t_ret;
	}
	;
};

/** Compose a pre_intersection_transform with a transform. */
template<typename T> class parallel_composition_operator<T, pre_intersection_transform> {
public:
	static continuous_set_transform::ptr implement(const T* t1,
			const pre_intersection_transform* t2) {
		assert(t1);
		assert(t2);
		return parallel_composition_operator<pre_intersection_transform, T>::implement(t2, t1);
	}
	;
};

/** Compose a pre_intersection_transform with a pre_intersection_transform. */
template<> class parallel_composition_operator<pre_intersection_transform,
		pre_intersection_transform> {
public:
	static continuous_set_transform::ptr implement(const pre_intersection_transform* t1,
			const pre_intersection_transform* t2) {
		assert(t1);
		assert(t2);
		return get_intersection(t1, t2);
	}
	;
};

/** Compose a post_intersection_transform with a transform. */
template<typename T> class parallel_composition_operator<post_intersection_transform, T> {
public:
	static continuous_set_transform::ptr implement(const post_intersection_transform* t1,
			const T* t2) {
		assert(t1);
		assert(t2);
		// t1 goes in back
		continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(new sequence_transform(
				t2, t1));
		return t_ret;
	}
	;
};

/** Compose a post_intersection_transform with a transform. */
template<typename T> class parallel_composition_operator<T, post_intersection_transform> {
public:
	static continuous_set_transform::ptr implement(const T* t1,
			const post_intersection_transform* t2) {
		assert(t1);
		assert(t2);
		return parallel_composition_operator<post_intersection_transform, T>::implement(t2, t1);
	}
	;
};

/** Compose a post_intersection_transform with a post_intersection_transform. */
template<> class parallel_composition_operator<post_intersection_transform,
		post_intersection_transform> {
public:
	static continuous_set_transform::ptr implement(const post_intersection_transform* t1,
			const post_intersection_transform* t2) {
		assert(t1);
		assert(t2);
		return get_intersection(t1, t2);
	}
	;
};

continuous_set_transform::ptr get_intersection(const pre_intersection_transform* t1,
		const pre_intersection_transform* t2) {
	assert(t1);
	assert(t2);
	// Intersect the two continuous_sets of the transforms
	continuous_set_ptr cset;
	cset = compute_intersection(t1->get_set(), t2->get_set());
	continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(
			new pre_intersection_transform(cset));
	assert(t_ret);
	return t_ret;
}

continuous_set_transform::ptr get_intersection(const post_intersection_transform* t1,
		const post_intersection_transform* t2) {
	assert(t1);
	assert(t2);
	// Intersect the two continuous_sets of the transforms
	continuous_set_ptr cset;
	cset = compute_intersection(t1->get_set(), t2->get_set());
	continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(
			new post_intersection_transform(cset));
	assert(t_ret);
	return t_ret;
}

continuous_set_transform::ptr get_intersection(const intersection_transform* t1,
		const intersection_transform* t2) {
	// Intersect the two continuous_sets of the transforms
	continuous_set_ptr cset;
	cset = compute_intersection(t1->get_set(), t2->get_set());
	continuous_set_transform::ptr t_ret = continuous_set_transform::ptr(
			new pre_intersection_transform(cset));
	assert(t_ret);
	return t_ret;
}

continuous_set_transform::ptr parallel_compose_sequence(const sequence_transform* pt1,
		const sequence_transform* pt2, bool try_swapping) {
	assert(pt1);
	assert(pt2);
	// test if the sequence is just one element
	continuous_set_transform::ptr t_ret = continuous_set_transform::ptr();
	unsigned int k1 = pt1->size();
	unsigned int k2 = pt2->size();
	if (k1 == 0) { // return a copy of pt2
		t_ret = continuous_set_transform::ptr(new sequence_transform(pt2->begin(), pt2->end()));
		return t_ret;
	} else if (k2 == 0) {// return a copy of pt1
		t_ret = continuous_set_transform::ptr(new sequence_transform(pt1->begin(), pt1->end()));
		return t_ret;
	} else if (k1 == 1 && k2 == 1) {
		return parallel_compose(pt1->front(), pt2->front());
	} else { // now either k1 or k2 have size >= 2
		// if they begin with a pre-intersection, merge them
		if (const pre_intersection_transform* pre_pt1
				= dynamic_cast<const pre_intersection_transform*>(pt1->begin()->get())) {
			if (const pre_intersection_transform* pre_pt2
					= dynamic_cast<const pre_intersection_transform*>(pt2->begin()->get())) {
				continuous_set_transform::ptr t_inter = get_intersection(pre_pt1, pre_pt2);
				// compose the rest of t1 with the rest of t2
				sequence_transform::ptr s1 = sequence_transform::ptr(new sequence_transform(
						++pt1->begin(), pt1->end()));
				sequence_transform::ptr s2 = sequence_transform::ptr(new sequence_transform(
						++pt2->begin(), pt2->end()));
				continuous_set_transform::ptr t_rest = parallel_composition_operator<
						sequence_transform, sequence_transform>::implement(s1.get(), s2.get());
				if (t_rest) {
					t_ret = continuous_set_transform::ptr(new sequence_transform(t_inter, t_rest));
					return t_ret;
				}
			}
		}
		// if they end with a post-intersection, merge them
		if (const post_intersection_transform* post_pt1
				= dynamic_cast<const post_intersection_transform*>(pt1->back().get())) {
			if (const post_intersection_transform* post_pt2
					= dynamic_cast<const post_intersection_transform*>(pt2->back().get())) {
				continuous_set_transform::ptr t_inter = get_intersection(post_pt1, post_pt2);
				// compose the rest of t1 with the rest of t2
				sequence_transform::ptr s1;
				sequence_transform::ptr s2;
				s1 = sequence_transform::ptr(new sequence_transform(pt1->begin(), --(pt1->end())));
				s2 = sequence_transform::ptr(new sequence_transform(pt2->begin(), --(pt2->end())));
				continuous_set_transform::ptr t_rest = parallel_composition_operator<
						sequence_transform, sequence_transform>::implement(s1.get(), s2.get());
				if (t_rest) {
					t_ret = continuous_set_transform::ptr(new sequence_transform(t_rest, t_inter));
					return t_ret;
				}
			}
		}
		// Otherwise try to sequentialize them one by one.
		// Get the variables used by t2, and find how many elements of t1 can be executed before modifying a variable that
		// is used by t2.
		variable_id_set used1, modif1, used2, modif2;
		pt1->front()->get_used_and_modif_variables(used1, modif1); // get the variables modified by the first transform of t1
		pt2->get_used_and_modif_variables(used2, modif2);
		if (set_is_disjoint<variable_id> (modif1, used2)) {
			// compose the rest of t1 with t2
			sequence_transform::ptr s1 = sequence_transform::ptr(new sequence_transform(
					++pt1->begin(), pt1->end()));
			continuous_set_transform::ptr t_rest = parallel_composition_operator<
					sequence_transform, sequence_transform>::implement(s1.get(), pt2);
			if (t_rest) {
				// return the first transform of t1 followed by the result of the composition
				t_ret = continuous_set_transform::ptr(new sequence_transform(pt1->front(), t_rest));
				return t_ret;
			}
		}
		// Last resort: try swapping t1 and t2
		if (try_swapping) {
			return parallel_compose_sequence(pt2, pt1, false); // false prevents from swapping again and again and again... 
		}
		return continuous_set_transform::ptr();
	}
}

}
