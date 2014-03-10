#include "core/continuous/continuous_set_transforms/sequence_transform.h"

/** Forward declarations */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
continuous_set_ptr compute_or_assign_relation_concatenation(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2);
}

#include "core/continuous/continuous_set.h"

namespace continuous {

/** Initialize to an empty sequence. */
sequence_transform::sequence_transform() {
}

/** Initialize with the singleton sequence (t1), flattening it if t1 is itself a sequence. */
sequence_transform::sequence_transform(const continuous_set_transform::const_ptr& t1) {
	push_back(t1);
}

/** Initialize with the sequence (t1,t2), joining them if one of them is a sequence. */
sequence_transform::sequence_transform(const continuous_set_transform::const_ptr& t1,
		const continuous_set_transform::const_ptr& t2) {
	push_back(t1);
	push_back(t2);
}

/** Initialize with the sequence (t1,t2,t3), joining them if one of them is a sequence. */
sequence_transform::sequence_transform(const continuous_set_transform::const_ptr& t1,
		const continuous_set_transform::const_ptr& t2,
		const continuous_set_transform::const_ptr& t3) {
	push_back(t1);
	push_back(t2);
	push_back(t3);
}

/** Initialize with the sequence (*it1,...,*it2). */
sequence_transform::sequence_transform(const const_iterator& it1, const const_iterator& it2) {
	for (const_iterator it = it1; it != it2; ++it) {
		transform_list.push_back(*it);
	}
}

sequence_transform::~sequence_transform() {
}

sequence_transform::const_iterator sequence_transform::begin() const {
	return transform_list.begin();
}

sequence_transform::const_iterator sequence_transform::end() const {
	return transform_list.end();
}

/** Returns the first element of the sequence. List must not be empty. */
continuous_set_transform::const_ptr sequence_transform::front() const {
	return transform_list.front();
}

/** Returns the last element of the sequence. List must not be empty. */
continuous_set_transform::const_ptr sequence_transform::back() const {
	return transform_list.back();
}

unsigned int sequence_transform::size() const {
	return transform_list.size();
}

void sequence_transform::push_front(const continuous_set_transform::const_ptr& t) {
	if (const sequence_transform* pt1 = dynamic_cast<const sequence_transform*>(t.get())) {
		for (sequence_transform::const_iterator it = pt1->begin(); it != pt1->end(); ++it) {
			transform_list.push_front(*it);
		}
	} else {
		transform_list.push_front(t);
	}
}

void sequence_transform::push_back(const continuous_set_transform::const_ptr& t) {
	if (const sequence_transform* pt1 = dynamic_cast<const sequence_transform*>(t.get())) {
		for (sequence_transform::const_iterator it = pt1->begin(); it != pt1->end(); ++it) {
			transform_list.push_back(*it);
		}
	} else {
		transform_list.push_back(t);
	}
}

relation_const_ptr sequence_transform::get_relation(continuous_set_const_ptr cset) const {
	if (transform_list.empty()) {
		return relation_const_ptr(cset->create_universe());
	} else {
		continuous_set_ptr p;
		for (transform_list_type::const_iterator it = transform_list.begin(); it
				!= transform_list.end(); ++it) {
			if (it == transform_list.begin()) {
				p = continuous_set_ptr((*it)->get_relation(cset)->clone());
			} else {
				p = compute_or_assign_relation_concatenation(p, (*it)->get_relation(cset));
			}
		}
		return p;
	}
}

void sequence_transform::get_used_and_modif_variables(variable_id_set& used_vars,
		variable_id_set& modif_vars) const {
	get_used_and_modif_variables(used_vars, modif_vars, transform_list.begin(),
			transform_list.end());
}

void sequence_transform::print(std::ostream& os) const {
	for (sequence_transform::const_iterator it = transform_list.begin(); it != transform_list.end(); ++it) {
		if (it != transform_list.begin())
			os << "; ";
		os << *it;
	}
}

/** \todo compute set_difference  more effiently without tmp_set. */
void sequence_transform::get_used_and_modif_variables(variable_id_set& used_vars,
		variable_id_set& modif_vars, const const_iterator& bit, const const_iterator& eit) const {
	// A variable is used in (bit,eit) if it is used in bit or used in (++bit,eit) and not modified in bit.
	// (It has to be used and not previously reassigned.)
	// A variable is modified in (bit,eit) if it is modified in bit or modified in (++bit,eit).
	// (It has to be const in all of them.)
	if (bit != eit) {
		// do the combination with the reset of them
		const_iterator incbit = bit;
		++incbit;
		if (incbit != eit) {
			get_used_and_modif_variables(used_vars, modif_vars, incbit, eit);
		}

		// get vars of bit
		variable_id_set bit_used_vars, bit_modif_vars;
		(*bit)->get_used_and_modif_variables(bit_used_vars, bit_modif_vars);

		// compute used_vars = bit_used_vars \cup (used_vars \setminus bit_modif_vars)
		variable_id_set tmp_set;
		set_difference(used_vars.begin(), used_vars.end(), bit_modif_vars.begin(),
				bit_modif_vars.end(), inserter(bit_used_vars, bit_used_vars.begin()));
		used_vars.swap(bit_used_vars);
		// compute modif_vars
		tmp_set = variable_id_set();
		set_union(bit_modif_vars.begin(), bit_modif_vars.end(), modif_vars.begin(),
				modif_vars.end(), inserter(tmp_set, tmp_set.begin()));
		modif_vars.swap(tmp_set); // replace the current set with the result computed above
	} else {
		used_vars = variable_id_set();
		modif_vars = variable_id_set();
	}
}

void sequence_transform::accept(const_visitor& d) const {
	d.dispatch(this);
}

}

