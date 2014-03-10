#include "core/symbolic_states/symbolic_state_collection.h"

#include "core/symbolic_states/symbolic_state.h"
#include "core/discrete/discrete_set.h"
#include "core/continuous/continuous_set.h"

namespace hybrid_automata {

/* need to declare static members */
symbolic_state_collection::output_format
		symbolic_state_collection::my_output_format;

symbolic_state_collection::ptr symbolic_state_collection::get_ptr() {
	symbolic_state_collection::ptr p = boost::enable_shared_from_this<
			symbolic_state_collection>::shared_from_this();
	return p;
}

symbolic_state_collection::const_ptr symbolic_state_collection::get_const_ptr() const {
	symbolic_state_collection::const_ptr p = boost::enable_shared_from_this<
			symbolic_state_collection>::shared_from_this();
	return p;
}

void symbolic_state_collection::copy(
		const symbolic_state_collection::const_ptr& sstate_set) {
	clear();
	for (const_iterator it = sstate_set->begin(); it != sstate_set->end(); ++it) {
		add((*it)->clone());
	};
}

variable_id_set symbolic_state_collection::get_variable_ids() const {
	variable_id_set vis;
	for (const_iterator it = begin(); it != end(); ++it) {
		variable_id_set vis2=(*it)->get_continuous_set()->get_variable_ids();
		std::set_union(vis.begin(), vis.end(), vis2.begin(), vis2.end(),
				inserter(vis, vis.begin()));
	}
	return vis;

}

bool symbolic_state_collection::is_disjoint_from(
		const symbolic_state_collection::ptr& sstate_set) const {
	/* Check if all elements of *this are disjoint from all elements of sstate_set. */
	for (const_iterator it = begin(); it != end(); ++it) {
		for (const_iterator jt = sstate_set->begin(); jt != sstate_set->end(); ++jt) {
			if (!(*it)->get_discrete_set()->is_disjoint_from(
					(*jt)->get_discrete_set())) {
				if (!(*it)->get_continuous_set()->is_disjoint_from(
						(*jt)->get_continuous_set()))
					return false;
			}
		}
	}
	return true;
}

void symbolic_state_collection::union_assign(
		const symbolic_state_collection::const_ptr& sstate_set) {
	for (const_iterator it = sstate_set->begin(); it != sstate_set->end(); ++it) {
		add(*it);
	}
}

void symbolic_state_collection::print(std::ostream& os) const {
	const symbolic_state_collection::output_format& of = get_output_format();
	os << of.preamble;
	if (begin() == end()) {
		os << of.empty_signal;
	} else {
		for (const_iterator it = begin(); it != end(); ++it) {
			if (it != begin())
				os << of.element_separator;
			os << *it;
		}
	}
	os << of.epilogue;
}

}
