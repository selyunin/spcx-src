#include "core/continuous/continuous_set_collection.h"

namespace continuous {

continuous_set_collection::output_format
		continuous_set_collection::my_output_format=
				continuous_set_collection::output_format();

continuous_set_collection::continuous_set_collection(continuous_set::ptr element) {
	my_container.push_back(element);
}

continuous_set_collection::~continuous_set_collection() {
}

continuous_set_collection* continuous_set_collection::clone() const {
	continuous_set_collection* new_coll = new continuous_set_collection();
	for (const_iterator it=begin(); it != end(); ++it) {
		continuous_set::ptr element_clone=continuous_set::ptr((*it)->clone());
		new_coll->push_back(element_clone);
	}
	return new_coll;
}

continuous_set_collection* continuous_set_collection::create_universe() const {
	continuous_set::ptr universe_element(default_element()->create_universe());
	return new continuous_set_collection(universe_element);
}

continuous_set_collection* continuous_set_collection::create_empty() const {
	continuous_set::ptr empty_element(default_element()->create_empty());
	return new continuous_set_collection(empty_element);
}

continuous_set::const_ptr continuous_set_collection::default_element() const {
	assert(size()>0);
	return *begin();
}

continuous_set_collection::iterator continuous_set_collection::begin() {
	return my_container.begin();
}

continuous_set_collection::iterator continuous_set_collection::end() {
	return my_container.end();
}

continuous_set_collection::const_iterator continuous_set_collection::begin() const {
	return my_container.begin();
}

continuous_set_collection::const_iterator continuous_set_collection::end() const {
	return my_container.end();
}

unsigned int continuous_set_collection::size() const {
	return my_container.size();
}

int continuous_set_collection::get_memory() const {
	int mem=0;
	for (const_iterator it=begin(); it != end(); ++it) {
		mem+=(*it)->get_memory();
	}
	return mem;
}

const variable_id_set& continuous_set_collection::get_variable_ids() const {
	static variable_id_set vis;
	vis = variable_id_set();
	for (const_iterator it=begin(); it != end(); ++it) {
		variable_id_set element_vars=(*it)->get_variable_ids();
		vis.insert(element_vars.begin(), element_vars.end());
	}
	return vis;
}

unsigned int continuous_set_collection::get_dim() const {
	return get_variable_ids().size();
}

math::tribool continuous_set_collection::is_empty() const {
	math::tribool result(true);
	for (const_iterator it=begin(); it != end() && math::maybe(result); ++it) {
		result = result && (*it)->is_empty();
	}
	return result;
}

void continuous_set_collection::embed_variables(const variable_id_set& id_set) {
	for (iterator it=begin(); it != end(); ++it) {
		(*it)->embed_variables(id_set);
	}
}

void continuous_set_collection::existentially_quantify_variables(
		const variable_id_set& id_set) {
	for (iterator it=begin(); it != end(); ++it) {
		(*it)->existentially_quantify_variables(id_set);
	}
}

math::tribool continuous_set_collection::element_wise_contains(
		continuous_set::const_ptr p) const {
	math::tribool result(false);
	for (const_iterator it=begin(); it != end() && math::maybe(!result); ++it) {
		result = result || (*it)->contains(p);
	}
	return result;
}

void continuous_set_collection::delete_if_contained_in(
		continuous_set::const_ptr p) {
	for (iterator it=begin(); it != end();) {
		if (math::definitely(p->contains(*it)))
			it=my_container.erase(it);
		else
			++it;
	}
}

/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
 */
void continuous_set_collection::reassign_primedness(unsigned int d, unsigned int p) {
	for (iterator it=begin(); it != end(); ++it) {
		(*it)->reassign_primedness(d, p);
	}
}

/** Increase the primedness of the variables with primedness of degree \p d by 1.
 * If d is 0, increase all. */
void continuous_set_collection::increase_primedness(unsigned int d) {
	for (iterator it=begin(); it != end(); ++it) {
		(*it)->increase_primedness(d);
	}
}

/** Decrease the primedness of the variables with primedness of degree \p d by 1.
 * If d is 0, decrease all. */
void continuous_set_collection::decrease_primedness(unsigned int d) {
	for (iterator it=begin(); it != end(); ++it) {
		(*it)->decrease_primedness(d);
	}
}

/** Insert p into the container. By default or with redundancy_check=true,
 * the element is only inserted if it is not contained already; also,
 * all elements contained in p are removed. */
void continuous_set_collection::insert(continuous_set::ptr p,
		bool redundancy_check) {
	if (my_container.empty() || !math::definitely(p->is_empty())) {
		/* test if p is already in the elements of *this.
		 The test is done element-wise since the actual operation involves
		 the difference operator, which is excessively costly. */
		if (!redundancy_check || !math::definitely(element_wise_contains(p))) {
			/* Test if p contains eny element of *this. */
			if (redundancy_check) {
				delete_if_contained_in(p);
			}
			push_back(p);
		}
	}
}

continuous_set_predicate::ptr continuous_set_collection::get_predicate() const {
	throw std::runtime_error("missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

void continuous_set_collection::accept(
		dispatching::dispatcher<continuous_set_typelist>& d) const {
	// The following is dangerous, since the dispatcher might be stateful.
	// Accepting it with more than one element could lead to wrong results.
	d.dispatch(this);
//	for (const_iterator it=begin(); it != end(); ++it) {
//		(*it)->accept(d);
//	}
}

void continuous_set_collection::print(std::ostream& os) const {
	output_format of=get_output_format();
	os << of.preamble;

	for (const_iterator it=begin(); it != end(); ++it) {
		if (it!=begin())
			os << of.element_separator;
		os << *it;
	}
	os << of.epilogue;
}

continuous_set_collection::continuous_set_collection() {
}

void continuous_set_collection::push_back(continuous_set::ptr p) {
	my_container.push_back(p);
}

}

